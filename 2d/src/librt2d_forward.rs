#![allow(dead_code)]
#![allow(unused)]

use std::{
    arch::x86_64::_MM_PERM_ABDA,
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
    rt_common::*,
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
pub struct ForwardRenderParams {
    pub width: i32,
    pub height: i32,
    pub stop_condition: StopCondition,
    pub recursion_limit: usize,
    pub lambda_samples: usize,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ForwardRenderer {
    render_params: ForwardRenderParams,
}

impl ForwardRenderer {
    pub fn new(render_params: ForwardRenderParams) -> Self {
        Self { render_params }
    }

    /// Compute list of ray bounces from origin
    /// return [(ray, hit), ... ]
    /// decrease `depth`
    /// until `depth` == 0
    pub fn trace_ray(
        &self,
        world: &World,
        ray: &Ray2d,
        lambda_idx: usize,
        lambda: f64,
        depth: usize,
        power: f32,
    ) -> Vec<(Ray2d, f32, Option<Hit2d>)> {
        if depth == 0 {
            return vec![];
        }

        let (obj, hit): (&Object, Hit2d) = {
            let first_hit: Option<(&Object, Hit2d)> = world
                .objects
                .iter()
                .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj, hit)))
                .min_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));
            let Some((obj, hit)) = first_hit else {
                return vec![(ray.clone(), power, None)];
            };
            (obj, hit.clone())
        };

        // delete this
        //return vec![(ray.clone(), power, Some(hit))];

        if hit.side == Side::Inside {
            match obj.mat {
                Material::Emissive { mut emission } => {
                    //emission.data[lambda_idx]
                    return vec![];
                }
                Material::Dielectric { ior, absorption } => {
                    let ior = ior.ior(lambda);
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(-hit.n, ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(-hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let new_power = power * f32::exp(-hit.t as f32 * absorption.data[lambda_idx]);
                    let next = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1, new_power);
                    next
                }
                //Material::SubSurfaceScattering { sigma_a, sigma_s } => {
                //    return vec![];
                //}
                _ => vec![],
            }
        } else {
            // outside
            match obj.mat {
                Material::Diffuse { absorption } => {
                    if absorption == Spectrum::default() {
                        // 100% absorption
                        return vec![];
                    }

                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::rand_hemisphere(p, hit.n);

                    let new_power = absorption.data[lambda_idx] * power;
                    if new_power > 0. {
                        let next =
                            self.trace_ray(world, &r, lambda_idx, lambda, depth - 1, new_power);
                        next
                    } else {
                        vec![]
                    }
                }
                Material::Emissive { mut emission } => {
                    //emission.data[lambda_idx]
                    return vec![];
                }
                Material::DirectionalEmissive {
                    emission: emission_color,
                    angle,
                    d,
                } => {
                    return vec![];
                    //if ray.dir.dot(-hit.n) >= angle {
                    //    let distance_coeff = 1. / (hit.t + d);
                    //    emission_color.data[lambda_idx] * distance_coeff as f32
                    //} else {
                    //    vec![]
                    //}
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, ray.dir.reflect(hit.n));

                    let new_power = power; // 100% reflective
                    let next = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1, new_power);
                    next
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

                    let new_power = power;
                    let next = self.trace_ray(world, &r, lambda_idx, lambda, depth - 1, power);
                    next
                }
                //Material::SubSurfaceScattering { absorption } => todo!(),
                _ => vec![],
            }
        }
    }

    /// renderer entry point
    pub fn global_render(&self, world: &World) -> RawImage {
        let mut raw_image = RawImage::new(self.render_params.width, self.render_params.height);

        // list emissive materials objects
        let emissive_objects: Vec<&Object> = world
            .objects
            .iter()
            .filter(|obj| match obj.mat {
                Material::Emissive { emission } => true,
                Material::DirectionalEmissive { emission, angle, d } => true,
                _ => false,
            })
            .collect();

        for obj in &emissive_objects {
            fn sample_emission(mat: Material) -> (usize, f64, f32) {
                let lambda_idx = 0;
                let lambda = 1.;
                let power = 1.;

                (lambda_idx, lambda, power)
            }

            let (lambda_idx, lambda, power) = sample_emission(obj.mat);

            for _ in 0..1000 {
                let ray = obj.shape.sample_out_ray();
                let bounce_list = self.trace_ray(world, &ray, lambda_idx, lambda, 5, power);

                for (ray, power, hit) in &bounce_list {
                    let t_max = hit.map_or(1000., |hit| hit.t);

                    use line_drawing::XiaolinWu;

                    let start = ray.origin;
                    let end = ray.origin + t_max * ray.dir;

                    for ((x, y), value) in
                        XiaolinWu::<f64, i32>::new((start.x, start.y), (end.x, end.y))
                    {
                        let pixel = IVec2::new(x, y);
                        let pixel_data = PixelData {
                            value: Vec3::ONE * value as f32,
                            weight: 1.,
                        };

                        let blending = Blending::Add;
                        raw_image.draw_pixel(pixel, pixel_data, blending);
                    }
                }
            }
        }

        // TODO
        // pick random object
        // pick random pos at edge of object
        // pick random dir
        // pixels list = trace ray()
        // color pixels on image
        // https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm

        raw_image
    }
}

impl Renderer for ForwardRenderer {
    fn endless_render(
        &self,
        world: World,
        tx: SyncSender<RenderProgress>,
        rx: Receiver<RenderCommand>,
    ) {
        println!(
            "Starting render. Stop condition={:?}",
            self.render_params.stop_condition
        );

        let mut merged_image = self.global_render(&world);
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
            let mut loop_image = self.global_render(&world);
            let global_render_elapsed = chrono.elapsed();
            println!("global_render: {:?}", global_render_elapsed);

            merged_image = merged_image + loop_image.clone();
            let elapsed = chrono.elapsed();
            println!("render loop = {:?}", elapsed);

            // find min/max for heatmap
            let mut heatmap = merged_image.map_pixel(|pixeldata| PixelData {
                weight: 1.,
                value: Vec3::splat(pixeldata.weight),
            });
            let max = heatmap
                .data
                .iter()
                .max_by(|a, b| a.value.x.total_cmp(&b.value.x))
                .unwrap()
                .value
                .x;
            let min = heatmap
                .data
                .iter()
                .min_by(|a, b| a.value.x.total_cmp(&b.value.x))
                .unwrap()
                .value
                .x;

            let colormap = colorgrad::preset::inferno();
            heatmap
                .map_pixel(|pixeldata| {
                    let t = (pixeldata.value.x - min) / max;
                    let color = colormap.at(t);
                    PixelData {
                        value: Vec3::new(color.r, color.g, color.b),
                        weight: 1.,
                    }
                })
                .convert_to_image(&ToneMappingMethod::Reinhard { param: 1. })
                .save("out/heatmap.png");

            //    .convert_to_image(&ToneMappingMethod::Reinhard);
            //heatmap.save("out/heatmap.png");

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
