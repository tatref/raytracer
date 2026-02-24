#![allow(dead_code)]
#![allow(unused)]

use std::{
    arch::x86_64::_MM_PERM_ABDA,
    collections::HashSet,
    f64::consts::PI,
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
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

pub fn annotate(raw_image: &mut RawImage, world: &World, p: DVec2) {
    let pixel_data = PixelData {
        value: Vec3::X,
        weight: 1.,
    };

    for i in 0..(world.render_params.width / 100) {
        for j in 0..5 {
            let pixel = IVec2::new(i * 100, j);
            let _ = raw_image.draw_pixel(pixel, pixel_data, Blending::Replace);
        }
    }

    for j in 0..(world.render_params.height / 100) {
        for i in 0..5 {
            let pixel = IVec2::new(i, j * 100);
            let _ = raw_image.draw_pixel(pixel, pixel_data, Blending::Replace);
        }
    }

    //let arcs = world.compute_arc_sections(p);

    //let r = 100.;
    //for i in 0..300 {
    //    let angle = i as f64 * 2. * PI / 300.;
    //    for arc in &arcs {
    //        if angle > arc.start && angle < arc.end {
    //            // inside
    //            let p2 = p - DVec2::from_angle(angle) * r;
    //            let px = p2.as_ivec2();

    //            raw_image.draw_pixel(px, Color::new(1., 0., 0.), Blending::Replace);
    //        }
    //    }
    //}
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Circle {
    pub center: DVec2,
    pub r: f64,
}

impl Circle {
    pub fn new(center: DVec2, r: f64) -> Self {
        Circle { center, r }
    }

    pub fn aabb(&self) -> Aabb {
        Aabb::new(self.center, DVec2::splat(self.r))
    }

    pub fn is_inside(&self, p: DVec2) -> bool {
        (self.center - p).length_squared() < self.r.powi(2)
    }

    pub fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
        let oc = ray.origin - self.center;
        let a = 1.0;
        let b = 2.0 * ray.dir.dot(oc);
        let c = oc.length_squared() - self.r.powi(2);
        let delta = b * b - 4.0 * a * c;

        if delta < 0.0 {
            return None;
        }

        let sqrt_delta = delta.sqrt();
        let t0 = (-b - sqrt_delta) / (2.0 * a);
        let t1 = (-b + sqrt_delta) / (2.0 * a);

        let t = if t0 > 0.0 {
            t0
        } else if t1 > 0.0 {
            t1
        } else {
            return None;
        };

        let p = ray.origin + t * ray.dir;
        let n = (p - self.center).normalize();
        let side = if self.is_inside(ray.origin) {
            Side::Inside
        } else {
            Side::Outside
        };

        Some(Hit2d { t, p, n, side })
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Bezier {
    control_points: Vec<DVec2>,
}

impl Bezier {
    pub fn new(control_points: &[DVec2]) -> Self {
        Self {
            control_points: control_points.iter().cloned().collect(),
        }
    }

    pub fn new_quadratic(a: DVec2, b: DVec2, c: DVec2) -> Self {
        let control_points = vec![a, b, c];
        Bezier::new(&control_points)
    }

    pub fn new_cubic(a: DVec2, b: DVec2, c: DVec2, d: DVec2) -> Self {
        let control_points = vec![a, b, c, d];
        Bezier::new(&control_points)
    }

    pub fn aabb(&self) -> Aabb {
        let mut max = DVec2::MAX;
        let mut min = DVec2::MIN;
        for p in &self.control_points {
            max = min.min(*p);
            min = max.min(*p);
        }

        let mid = min.midpoint(max);
        let half_size = (max - min) / 2.;

        Aabb::new(mid, half_size)
    }

    pub fn at(&self, t: f64) -> DVec2 {
        let order = self.control_points.len();

        let mut previous_segments = self.control_points.clone();

        for _ in 1..order {
            let mut new_segments = Vec::new();

            for tuple in previous_segments.windows(2) {
                let p = tuple[0].lerp(tuple[1], t);
                new_segments.push(p);
            }

            previous_segments = new_segments;
        }

        assert_eq!(previous_segments.len(), 1);
        previous_segments[0]
    }

    pub fn as_segments(&self, resolution: usize) -> Vec<Segment> {
        let segments: Vec<Segment> = (0..resolution)
            .map(|t| t as f64 / resolution as f64)
            .tuple_windows()
            .map(|(t0, t1)| {
                let p0 = self.at(t0);
                let p1 = self.at(t1);

                Segment::new(p0, p1)
            })
            .collect();

        segments
    }
}

#[derive(Debug)]
pub struct ArcSection {
    pub start: f64,
    pub end: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Segment {
    a: DVec2,
    b: DVec2,
    n: DVec2,
}

impl Segment {
    pub fn new(a: DVec2, b: DVec2) -> Self {
        let n = (b - a).perp().normalize();
        Segment { a, b, n }
    }

    pub fn flip(self) -> Self {
        Segment::new(self.b, self.a)
    }

    pub fn aabb(&self) -> Aabb {
        let max = self.a.max(self.b);
        let min = self.a.min(self.b);

        let mid = min.midpoint(max);
        let half_size = (max - min) / 2.;

        Aabb::new(mid, half_size)
    }

    // https://en.wikipedia.org/wiki/Circular-arc_graph
    pub fn get_arc(&self, p: DVec2) -> ArcSection {
        let p_a = self.a - p;
        let p_b = self.b - p;

        let angle_a = p_a.to_angle() + PI;
        let angle_b = p_b.to_angle() + PI;

        ArcSection {
            start: angle_a,
            end: angle_b,
        }
    }

    pub fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
        let r = ray.dir;
        let s = self.b - self.a;

        let denom = r.perp_dot(s);
        // Si denom == 0 → parallèle (pas d’intersection)
        if denom.abs() < f64::EPSILON {
            return None;
        }

        let diff = self.a - ray.origin;

        let t = diff.perp_dot(s) / denom;
        let u = diff.perp_dot(r) / denom;

        // Conditions :
        // - t >= 0 : intersection en avant du rayon
        // - 0 <= u <= 1 : intersection sur le segment
        if t >= 0.0 && u >= 0.0 && u <= 1.0 {
            let p = ray.origin + t * r;
            let side = if self.n.dot(ray.dir) < 0.0 {
                Side::Outside
            } else {
                Side::Inside
            };
            Some(Hit2d {
                t,
                p,
                n: self.n,
                side,
            })
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Strip {
    strip: Vec<Segment>,
}

impl Strip {
    pub fn new(strip: &[Segment]) -> Self {
        Self {
            strip: strip.iter().cloned().collect(),
        }
    }

    pub fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
        let mut hits: Vec<Hit2d> = self
            .strip
            .iter()
            .filter_map(|segment| segment.hit(ray))
            .collect();

        hits.sort_by(|a, b| a.t.total_cmp(&b.t));

        hits.first().cloned()
    }

    pub fn aabb(&self) -> Aabb {
        let mut max = DVec2::NEG_INFINITY;
        let mut min = DVec2::INFINITY;

        for seg in &self.strip {
            max = max.max(seg.a.max(seg.b));
            min = min.min(seg.a.min(seg.b));
        }

        let mid = min.midpoint(max);
        let half_size = (max - min) / 2.;

        Aabb::new(mid, half_size)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Shape {
    Circle(Circle),
    Segment(Segment),
    Strip(Strip),
    //Bezier(Bezier),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Inside,
    Outside,
}

#[derive(Copy, Clone, Debug)]
pub struct Hit2d {
    pub(crate) t: f64,
    pub(crate) p: DVec2,
    pub(crate) n: DVec2,
    pub(crate) side: Side,
}

pub struct Ray2d {
    pub origin: DVec2,
    pub dir: DVec2,
}

impl Ray2d {
    pub fn new(origin: DVec2, dir: DVec2) -> Self {
        Self { origin, dir }
    }

    pub fn rand(origin: DVec2) -> Self {
        let angle: f64 = rand::random();
        let dir = DVec2::from_angle(2. * PI * angle);

        Self { origin, dir }
    }

    pub fn rand_hemisphere(origin: DVec2, n: DVec2) -> Self {
        let dir = loop {
            let v = DVec2::new(
                rand::random::<f64>() * 2.0 - 1.0,
                rand::random::<f64>() * 2.0 - 1.0,
            );
            if v.length_squared() <= 1.0 {
                let v = v.normalize();
                if v.dot(n) > 0.0 {
                    break v;
                }
            }
        };

        Self { origin, dir }
    }

    pub fn hit(&self, shape: &Shape) -> Option<Hit2d> {
        match shape {
            Shape::Circle(circle) => circle.hit(self),
            Shape::Segment(segment) => segment.hit(self),
            Shape::Strip(strip) => strip.hit(self),
            //Shape::Bezier(bezier) => bezier.hit(self),
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum Material {
    Emissive {
        /// emission color
        emission: Spectrum,
    },
    DirectionalEmissive {
        /// emission color
        emission: Spectrum,
        angle: f64,
        /// inner radius
        /// for better fall off
        d: f64,
    },
    Diffuse {
        absorption: Spectrum,
    },
    Reflective,
    Dielectric {
        ior: Ior,
    },
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum Ior {
    Simple(f64),
    /// https://en.wikipedia.org/wiki/Cauchy%27s_equation
    Cauchy {
        a: f64,
        b: f64,
    },
}
impl Ior {
    pub fn ior(&self, lambda: f64) -> f64 {
        let lambda = lambda / 1000.;
        match self {
            Ior::Simple(ior) => *ior,
            Ior::Cauchy { a, b } => *a + *b / lambda.powi(2),
        }
    }
}

impl Material {
    pub fn emissive(spectrum: Spectrum) -> Self {
        let color_at_surface = spectrum;
        Material::Emissive {
            emission: color_at_surface,
        }
    }
    pub fn directional_emissive_at(d: f64, angle: f64, emission: Spectrum) -> Self {
        let k = emission * d as f32;
        Material::DirectionalEmissive {
            emission: k,
            d,
            angle,
        }
    }

    pub fn diffuse(absorption: Spectrum) -> Self {
        Material::Diffuse { absorption }
    }

    pub fn dieletric(ior: f64) -> Self {
        Material::Dielectric {
            ior: Ior::Simple(ior),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Object {
    pub shape: Shape,
    pub mat: Material,
}

impl Object {
    pub fn new(shape: Shape, mat: Material) -> Self {
        Self { shape, mat }
    }

    pub fn segment(a: DVec2, b: DVec2, mat: Material) -> Self {
        let shape = Shape::Segment(Segment::new(a, b));
        Self { shape, mat }
    }

    pub fn bezier(points: &[DVec2], resolution: usize, mat: Material) -> Self {
        let strip = Bezier::new(points).as_segments(resolution);
        let strip = Strip { strip };

        Object {
            shape: Shape::Strip(strip),
            mat,
        }
    }

    pub fn from_points(points: &[DVec2], mat: Material) -> Self {
        assert!(points.len() >= 2);

        let mut strip: Vec<Segment> = points
            .windows(2)
            .map(|points| Segment::new(points[0], points[1]))
            .collect();
        let closing = Segment::new(points[points.len() - 1], points[0]);
        strip.push(closing);

        Object {
            shape: Shape::Strip(Strip { strip }),
            mat,
        }
    }
    pub fn aabb(&self) -> Aabb {
        match &self.shape {
            Shape::Circle(c) => c.aabb(),
            Shape::Segment(s) => s.aabb(),
            Shape::Strip(s) => s.aabb(),
        }
    }

    pub fn draw_aabb(&self, raw_image: &mut RawImage, color: Color) {
        let aabb = self.aabb();
        aabb.draw(raw_image, color);
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct RenderParams {
    pub width: i32,
    pub height: i32,
    pub spp: usize,
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

pub struct RenderProgress {
    pub loops: usize,
    pub raw_image: RawImage,
}

#[derive(Debug)]
pub enum RenderCommand {
    StopRender,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct World {
    pub render_params: RenderParams,
    //light: Light,
    pub objects: Vec<Object>,
    quadtree: Node,
}

impl World {
    pub fn new(objects: Vec<Object>, render_params: RenderParams) -> Self {
        let mut quadtree = Node::new(Aabb::new(
            DVec2::new(800. / 2., 600. / 2.),
            DVec2::new(800. / 2., 600. / 2.),
        ));

        Self {
            objects,
            quadtree,
            render_params,
        }
    }

    fn save(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::File::create(path)?;
        rmp_serde::encode::write(&mut f, self)?;

        Ok(())
    }

    pub fn compute_arc_sections(&self, p: DVec2) -> Vec<ArcSection> {
        let mut arcs = Vec::new();

        for obj in &self.objects {
            match &obj.shape {
                Shape::Segment(segment) => arcs.push(segment.get_arc(p)),
                _ => (),
            }
        }

        arcs
    }

    pub fn compute_pixel(&self, pixel: IVec2, spp: usize, recursion_limit: usize) -> Spectrum {
        let mut total_spectrum = Spectrum::default();

        for _ in 0..spp {
            let p = pixel.as_dvec2() + rand::random::<DVec2>() - DVec2::ONE / 2.;
            let ray = Ray2d::rand(p);

            let mut spectrum_spp = Spectrum::default();
            for _ in 0..self.render_params.lambda_samples {
                let (lambda_idx, lambda) = Spectrum::rand_lambda();
                let spectrum_lambda_sample =
                    self.trace_ray(&ray, lambda_idx, lambda, recursion_limit);
                spectrum_spp += spectrum_lambda_sample;
            }
            spectrum_spp =
                spectrum_spp / self.render_params.lambda_samples as f32 * SPECTRUM_SAMPLES as f32;

            total_spectrum += spectrum_spp;
        }
        total_spectrum / spp as f32
    }

    pub fn trace_ray(&self, ray: &Ray2d, lambda_idx: usize, lambda: f64, depth: usize) -> Spectrum {
        if depth == 0 {
            return Spectrum::default();
        }

        let (obj, hit): (&Object, Hit2d) = if self.render_params.use_quadtree {
            // Quadtree version
            let hit: Option<(&Object, Hit2d)> = self.quadtree.hit(ray);
            let Some((obj, hit)) = hit else {
                return Spectrum::default();
            };
            (obj, hit)
        } else {
            // Vec version
            let mut hits: Vec<(&Object, Hit2d)> = self
                .objects
                .iter()
                .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj, hit)))
                .collect();
            hits.sort_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));
            let first_hit = hits.first();
            let Some((obj, hit)) = first_hit else {
                return Spectrum::default();
            };
            (obj, hit.clone())
        };

        if hit.side == Side::Inside {
            match obj.mat {
                Material::Emissive { mut emission } => {
                    for idx in 0..SPECTRUM_SAMPLES {
                        if idx != lambda_idx {
                            emission.data[idx] = 0.;
                        }
                    }
                    emission
                }
                Material::Dielectric { ior } => {
                    let ior = ior.ior(lambda);
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(-hit.n, ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(-hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let spectrum = self.trace_ray(&r, lambda_idx, lambda, depth - 1);
                    spectrum
                }
                _ => Spectrum::default(),
            }
        } else {
            // outside
            match obj.mat {
                Material::Diffuse { absorption } => {
                    if absorption == Spectrum::default() {
                        // 100% absorption
                        return Spectrum::default();
                    }

                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::rand_hemisphere(p, hit.n);

                    let light = self.trace_ray(&r, lambda_idx, lambda, depth - 1);

                    absorption * light
                }
                Material::Emissive { mut emission } => {
                    for idx in 0..SPECTRUM_SAMPLES {
                        if idx != lambda_idx {
                            emission.data[idx] = 0.;
                        }
                    }
                    emission
                }
                Material::DirectionalEmissive {
                    emission: emission_color,
                    angle,
                    d,
                } => {
                    if ray.dir.dot(-hit.n) >= angle {
                        let distance_coeff = 1. / (hit.t + d);
                        emission_color * distance_coeff as f32
                    } else {
                        Spectrum::default()
                    }
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, ray.dir.reflect(hit.n));

                    let col = self.trace_ray(&r, lambda_idx, lambda, depth - 1);
                    col
                }
                Material::Dielectric { ior } => {
                    let ior = ior.ior(lambda);
                    let p = hit.p - hit.n * 100000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(hit.n, 1. / ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let spectrum = self.trace_ray(&r, lambda_idx, lambda, depth - 1);
                    spectrum
                }
            }
        }
    }

    pub fn render_column(
        &self,
        render_params: &RenderParams,
        i: i32,
    ) -> Vec<(i32, (Spectrum, usize))> {
        let pixels: Vec<(i32, (Spectrum, usize))> = (0..render_params.height)
            .into_par_iter()
            .map(|j| {
                let pixel = IVec2::new(i, j);
                let spectrum =
                    self.compute_pixel(pixel, render_params.spp, render_params.recursion_limit);
                (j, (spectrum, render_params.spp))
            })
            .collect();
        pixels
    }

    // single uniform pass for all pixels
    pub fn global_render(&self) -> RawImage {
        let mut raw_image = RawImage::new(self.render_params.width, self.render_params.height);

        // parallel v1
        for i in 0..self.render_params.width {
            let pixels = self.render_column(&self.render_params, i);

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

    pub fn build_quadtree(&mut self) {
        if self.render_params.use_quadtree {
            let chrono = std::time::Instant::now();
            for obj in self.objects.iter().cloned() {
                self.quadtree.insert(obj);
            }
            let elapsed = chrono.elapsed();
            println!("build quadtree time: {:?}", elapsed);
        }
    }

    pub fn endless_render(&self, tx: SyncSender<RenderProgress>, rx: Receiver<RenderCommand>) {
        let mut merged_image = self.global_render();
        let render_progress = RenderProgress {
            loops: 0,
            raw_image: merged_image.clone(),
        };
        if let Err(e) = tx.try_send(render_progress) {
            println!("Can't send render_progress: {:?}", e);
        }

        let mut i = 1;
        let total_chrono = std::time::Instant::now();
        loop {
            println!("loop {}", i);

            let render_progress = RenderProgress {
                loops: i,
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
            let mut loop_image = self.global_render();
            let global_render_elapsed = chrono.elapsed();
            println!("global_render: {:?}", global_render_elapsed);

            merged_image = merged_image + loop_image.clone();
            match &self.render_params.denoiser {
                Some(denoiser) => {
                    for pass in 0..denoiser.passes {
                        let chrono = std::time::Instant::now();
                        self.denoise(&self.render_params, &mut merged_image, &loop_image, pass);
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
                .convert_to_image(&ToneMappingMethod::Reinhard)
                .save("out/heatmap.png");

            //    .convert_to_image(&ToneMappingMethod::Reinhard);
            //heatmap.save("out/heatmap.png");

            let chrono = std::time::Instant::now();

            i += 1;
        }
        let total_render_time = total_chrono.elapsed();
        println!("total render time: {:?}", total_render_time);
    }

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
        render_params: &RenderParams,
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
            .convert_to_image(&ToneMappingMethod::Reinhard)
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
            .convert_to_image(&ToneMappingMethod::Reinhard)
            .save(format!("out/mask_{}.png", pass))
        {
            println!("Can't save mask: {:?}", e);
        }

        let xxx: Vec<_> = recompute_pixels
            .par_iter()
            .map(|&pixel| {
                let color = self.compute_pixel(
                    pixel,
                    (render_params.spp as f32 * denoiser.oversampling_factor) as usize,
                    render_params.recursion_limit,
                );
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
