#![allow(dead_code)]
#![allow(unused)]

use std::{
    collections::HashSet,
    f64::consts::PI,
    fmt::Display,
    io::Write,
    os::raw,
    sync::mpsc::{Receiver, Sender, SyncSender},
};

use crate::{
    Color,
    aabb::{Aabb, Node},
    img::{Blending, PixelData, RawImage, ToneMappingMethod},
    spectrum::{SPECTRUM_SAMPLES, Spectrum},
};
use colorgrad::Gradient;
use glam::{DVec2, IVec2, Vec3};
use itertools::{Itertools, iproduct};
use rand_distr::Distribution;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use enum2egui::{Gui, GuiInspect};
use enum2str::EnumStr;

use crate::rt_common::*;

#[derive(Serialize, Deserialize, Clone, Copy, Debug, Gui)]
pub enum StopCondition {
    Endless,
    MaxLoops(usize),
}

impl Display for StopCondition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ReverseRenderParams {
    pub width: i32,
    pub height: i32,
    pub spp: usize,
    pub stop_condition: StopCondition,
    pub recursion_limit: usize,
    pub lambda_samples: usize,
    pub denoiser: Option<Denoiser>,
    pub use_quadtree: bool,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Denoiser {
    pub top: f32,
    pub mask_size: i32,
    pub oversampling_factor: f32,
    pub passes: usize,
}

pub struct ReverseRenderer {
    render_params: ReverseRenderParams,
}

impl ReverseRenderer {
    pub fn new(render_params: &ReverseRenderParams) -> Self {
        Self {
            render_params: render_params.clone(),
        }
    }
    //fn save(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    //    let mut f = std::fs::File::create(path)?;
    //    rmp_serde::encode::write(&mut f, self)?;

    //    Ok(())
    //}

    pub fn compute_pixel(
        &self,
        world: &World,
        pixel: IVec2,
        render_params: &ReverseRenderParams,
    ) -> Spectrum {
        let mut total_spectrum = Spectrum::default();

        for _ in 0..render_params.spp {
            let p = pixel.as_dvec2() + rand::random::<DVec2>() - DVec2::ONE / 2.;
            let ray = Ray2d::rand(p);

            let mut spectrum_spp = Spectrum::default();
            for _ in 0..render_params.lambda_samples {
                let (lambda_idx, lambda) = Spectrum::rand_lambda();
                let spectrum_lambda_sample = self.trace_ray(
                    world,
                    &ray,
                    lambda_idx,
                    lambda,
                    render_params.recursion_limit,
                );
                spectrum_spp.data[lambda_idx] += spectrum_lambda_sample;
            }
            spectrum_spp =
                spectrum_spp / render_params.lambda_samples as f32 * SPECTRUM_SAMPLES as f32;

            total_spectrum += spectrum_spp;
        }
        total_spectrum / render_params.spp as f32
    }

    pub fn trace_ray(
        &self,
        world: &World,
        ray: &Ray2d,
        lambda_idx: usize,
        lambda: f64,
        depth: usize,
    ) -> f32 {
        if depth == 0 {
            return f32::default();
        }

        let use_quadtree = false;

        let (obj, hit): (&Object, Hit2d) = if use_quadtree {
            // Quadtree version
            let hit: Option<(&Object, Hit2d)> = world.quadtree.hit(ray);
            let Some((obj, hit)) = hit else {
                return f32::default();
            };
            (obj, hit)
        } else {
            // Vec version
            let first_hit: Option<(&Object, Hit2d)> = world
                .objects
                .iter()
                .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj, hit)))
                .min_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));
            let Some((obj, hit)) = first_hit else {
                return f32::default();
            };
            (obj, hit.clone())
        };

        if hit.side == Side::Inside {
            match obj.mat {
                Material::Emissive { mut emission } => emission.data[lambda_idx],
                Material::Dielectric { ior, absorption } => {
                    let ior = ior.ior(lambda);
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(-hit.n, ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(-hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let power = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1);
                    power * f32::exp(-hit.t as f32 * absorption.data[lambda_idx])
                }
                Material::SubSurfaceScattering { sigma_a, sigma_s } => {
                    return f32::default();

                    //TODO

                    let sigma_t = sigma_a.data[lambda_idx] + sigma_s.data[lambda_idx];
                    if sigma_t == 0. {
                        return f32::default();
                    }

                    let bounces = 10;
                    let mut absorb = 1.;
                    let mut p;
                    let mut r;

                    while bounces > 0 {
                        p = hit.p + hit.n * 10000. * f64::EPSILON;
                        r = Ray2d::rand(p);
                        fn find_next_hit_distance(r: &Ray2d, world: &World) -> f64 {
                            0.
                        }
                        let d_surface = find_next_hit_distance(&r, world);

                        let d2: f64 = rand_distr::Exp::new(sigma_t as f64)
                            .unwrap()
                            .sample(&mut rand::rng());

                        if d_surface < d2 {
                            // ouside
                            todo!();
                            break;
                        } else {
                            // still inside, continue random walk
                            absorb *= sigma_s.data[lambda_idx] / sigma_t;
                        }
                        bounces -= 1;
                    }
                    if bounces == 0 {
                        // black
                        return f32::default();
                    }

                    let power = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1);
                    return power * absorb;
                }
                _ => f32::default(),
            }
        } else {
            // outside
            match obj.mat {
                Material::Diffuse { absorption } => {
                    if absorption == Spectrum::default() {
                        // 100% absorption
                        return f32::default();
                    }

                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::rand_hemisphere(p, hit.n);

                    let light = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1);

                    absorption.data[lambda_idx] * light
                }
                Material::Emissive { mut emission } => emission.data[lambda_idx],
                Material::DirectionalEmissive {
                    emission: emission_color,
                    angle,
                    d,
                } => {
                    if ray.dir.dot(-hit.n) >= angle {
                        let distance_coeff = 1. / (hit.t + d);
                        emission_color.data[lambda_idx] * distance_coeff as f32
                    } else {
                        f32::default()
                    }
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, ray.dir.reflect(hit.n));

                    let col = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1);
                    col
                }
                Material::Dielectric { ior, absorption } => {
                    let ior = ior.ior(lambda);
                    let p = hit.p - hit.n * 100000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(hit.n, 1. / ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let power = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1);
                    power
                }
                //Material::SubSurfaceScattering { absorption } => todo!(),
                _ => f32::default(),
            }
        }
    }

    pub fn render_column(
        &self,
        world: &World,
        render_params: &ReverseRenderParams,
        i: i32,
    ) -> Vec<(i32, (Spectrum, usize))> {
        let pixels: Vec<(i32, (Spectrum, usize))> = (0..render_params.height)
            .into_par_iter()
            .map(|j| {
                let pixel = IVec2::new(i, j);
                let spectrum = self.compute_pixel(world, pixel, &render_params);
                (j, (spectrum, render_params.spp))
            })
            .collect();
        pixels
    }

    // single uniform pass for all pixels
    pub fn global_render(&self, world: &World, render_params: &ReverseRenderParams) -> RawImage {
        let mut raw_image = RawImage::new(render_params.width, render_params.height);

        // parallel v1
        for i in 0..render_params.width {
            let pixels = self.render_column(world, render_params, i);

            for (j, (spectrum, weight)) in pixels {
                let pixel = IVec2::new(i, j);
                let pixel_data = PixelData {
                    value: spectrum.to_dvec3(),
                    weight: weight as f32,
                };
                raw_image
                    .draw_pixel(pixel, pixel_data, Blending::Replace)
                    .unwrap();
            }
        }

        raw_image
    }

    //pub fn render(&self) -> RawImage {
    //    //let aabb_color = Color::new(1., 0., 1.);
    //    //self.quadtree.draw(&mut raw_image, aabb_color);

    //    let mut raw_image = self.global_render();

    //    let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    //    image.save("raw.png").unwrap();

    //    match &self.render_params.denoiser {
    //        Some(denoiser) => {
    //            for pass in 0..denoiser.passes {
    //                let chrono = std::time::Instant::now();
    //                self.denoise(&self.render_params, &mut raw_image, pass);
    //                let elapsed = chrono.elapsed();
    //                println!("denoise pass={}/{}: {:?}", pass, denoiser.passes, elapsed);
    //                //let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    //                //image.save(format!("out/raw_{}.png", pass)).unwrap();
    //            }
    //            raw_image
    //        }
    //        None => raw_image,
    //    }
    //}

    /// denoise algorithm:
    /// Take 2 images as input: A with a lot of samples (reference image), B with less samples
    /// 1) calculate difference between the images, this highlights the noisy areas
    /// 2) perform a laplace transform
    /// 3) sort the pixels with the highest laplacian values
    /// 4) add a mask surrounding the top pixels
    /// 5) merge the masks together
    /// 6) recompute the masked pixels
    fn denoise(
        &self,
        world: &World,
        render_params: &ReverseRenderParams,
        merged_image: &mut RawImage,
        loop_image: &RawImage,
        pass: usize,
    ) {
        let denoiser = render_params.denoiser.as_ref().unwrap();
        let diff_image = (merged_image.clone() - loop_image.clone()).abs();

        let laplace_kernel = [[-1., -1., -1.], [-1., 8., -1.], [-1., -1., -1.]];
        let laplace_image = diff_image
            .convolution(laplace_kernel)
            .map_pixel(|p| PixelData {
                value: Vec3::splat(p.value.length()),
                weight: 1.,
            });
        if let Err(e) = laplace_image
            .convert_to_image(&ToneMappingMethod::Reinhard { param: 1. })
            .save(format!("out/laplace_{}.png", pass))
        {
            println!("can't save laplace: {:?}", e);
        }

        let mut sorted_pixels: Vec<(usize, &PixelData)> =
            laplace_image.data.iter().enumerate().collect();
        sorted_pixels.sort_by(|a, b| a.1.value.length().total_cmp(&b.1.value.length()));
        sorted_pixels.reverse();

        let total_pixels = render_params.width * render_params.height;
        let top = (total_pixels as f32 * denoiser.top) as usize;
        let recompute_indexes: Vec<usize> = sorted_pixels[0..top]
            .iter()
            .map(|(idx, _val)| *idx)
            .collect();

        let mut recompute_mask = merged_image.clone();
        let recompute_pixels: HashSet<IVec2> = recompute_indexes
            .iter()
            .map(|idx| {
                let mut pixels = Vec::new();
                let pixel = recompute_mask.idx_to_pixel(*idx).unwrap();
                let mask_size = denoiser.mask_size;
                for (i, j) in iproduct!(-mask_size..=mask_size, -mask_size..=mask_size) {
                    pixels.push(pixel + IVec2::new(i, j));
                }

                pixels
            })
            .flatten()
            .collect();

        for pixel in recompute_pixels.iter() {
            let _ = recompute_mask.draw_pixel(
                *pixel,
                PixelData {
                    weight: 1.,
                    value: Vec3::new(1., 0., 1.),
                },
                Blending::Replace,
            );
        }
        if let Err(e) = recompute_mask
            .convert_to_image(&ToneMappingMethod::Reinhard { param: 1. })
            .save(format!("out/mask_{}.png", pass))
        {
            println!("Can't save mask: {:?}", e);
        }

        let xxx: Vec<_> = recompute_pixels
            .par_iter()
            .map(|&pixel| {
                let color = self.compute_pixel(world, pixel, render_params);
                (pixel, color)
            })
            .collect();
        for (pixel, spectrum) in xxx {
            let pixel_data = PixelData {
                value: spectrum.to_dvec3(),
                weight: (render_params.spp as f32 * denoiser.oversampling_factor),
            };
            let _ = merged_image.draw_pixel(pixel, pixel_data, Blending::Add);
        }
    }
}

impl Renderer for ReverseRenderer {
    fn endless_render(
        &self,
        world: &World,
        tx: SyncSender<RenderProgress>,
        rx: Receiver<RenderCommand>,
    ) {
        println!(
            "Starting render. Stop condition={:?}",
            self.render_params.stop_condition
        );

        let mut merged_image = self.global_render(world, &self.render_params);
        let render_progress = RenderProgress {
            loops: 0,
            raw_image: merged_image.clone(),
        };
        if let Err(e) = tx.try_send(render_progress) {
            println!("Can't send render_progress: {:?}", e);
        }

        let mut loop_n = 1;
        let total_chrono = std::time::Instant::now();
        loop {
            println!("loop {}", loop_n);

            let render_progress = RenderProgress {
                loops: loop_n,
                raw_image: merged_image.clone(),
            };
            if let Err(e) = tx.try_send(render_progress) {
                println!("Can't send render_progress: {:?}", e);
            }
            if let Ok(render_command) = rx.try_recv() {
                println!("Received RenderCommand: {:?}", render_command);
                break;
            }

            let chrono = std::time::Instant::now();
            let mut loop_image = self.global_render(world, &self.render_params);
            let global_render_elapsed = chrono.elapsed();
            println!("global_render: {:?}", global_render_elapsed);

            merged_image = merged_image + loop_image.clone();
            match &self.render_params.denoiser {
                Some(denoiser) => {
                    for pass in 0..denoiser.passes {
                        let chrono = std::time::Instant::now();
                        self.denoise(
                            world,
                            &self.render_params,
                            &mut merged_image,
                            &loop_image,
                            pass,
                        );
                        println!(
                            "denoise {}/{}: {:?}",
                            pass,
                            denoiser.passes,
                            chrono.elapsed()
                        );
                    }
                }
                None => (),
            }
            let elapsed = chrono.elapsed();
            println!("render loop = {:?}", elapsed);

            merged_image.to_heatmap().save("out/heatmap.png");

            let chrono = std::time::Instant::now();

            loop_n += 1;
            match self.render_params.stop_condition {
                StopCondition::Endless => (),
                StopCondition::MaxLoops(max_loop) => {
                    if loop_n > max_loop {
                        break;
                    }
                }
            }
        }
        let total_render_time = total_chrono.elapsed();
        println!("Finished!");
        println!("total render time: {:?}", total_render_time);
    }
}
