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

pub fn annotate(raw_image: &mut RawImage, camera: &Camera, p: DVec2) {
    let pixel_data = PixelData {
        value: Vec3::X,
        weight: 1.,
    };

    for i in 0..(camera.size.x as i32 / 100) {
        for j in 0..5 {
            let pixel = IVec2::new(i * 100, j);
            let _ = raw_image.draw_pixel(pixel, pixel_data, Blending::Replace);
        }
    }

    for j in 0..(camera.size.y as i32 / 100) {
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

    pub fn sample_out_ray(&self) -> Ray2d {
        let angle = rand::random::<f64>() * 2. * PI;
        let n = DVec2::from_angle(angle);
        let delta = 0.0001; // TODO: proper fix
        let p = self.center + (self.r + delta) * n;

        Ray2d::new(p, n)
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

    pub fn sample_out_ray(&self) -> Ray2d {
        todo!()
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

    pub fn sample_out_ray(&self) -> Ray2d {
        todo!()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Shape {
    Circle(Circle),
    Segment(Segment),
    Strip(Strip),
    //Bezier(Bezier),
}
impl Shape {
    pub fn sample_out_ray(&self) -> Ray2d {
        match self {
            Shape::Circle(inner) => inner.sample_out_ray(),
            Shape::Segment(inner) => inner.sample_out_ray(),
            Shape::Strip(inner) => inner.sample_out_ray(),
        }
    }
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

#[derive(Clone)]
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
        /// white: no absorption, black: total absorption
        /// https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law
        absorption: Spectrum,
    },
    //SubSurfaceScattering {
    //    sigma_a: Spectrum,
    //    sigma_s: Spectrum,
    //},
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

    pub fn dieletric(ior: f64, absorption: Spectrum) -> Self {
        Material::Dielectric {
            ior: Ior::Simple(ior),
            absorption,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Camera {
    pub center: DVec2,
    pub size: DVec2,
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

#[derive(Clone, Serialize, Deserialize)]
pub struct World {
    pub camera: Camera,
    pub objects: Vec<Object>,
    pub quadtree: Node,
}

impl World {
    pub fn new(objects: Vec<Object>, camera: &Camera) -> Self {
        let mut quadtree = Node::new(Aabb::new(
            DVec2::new(800. / 2., 600. / 2.),
            DVec2::new(800. / 2., 600. / 2.),
        ));

        Self {
            objects,
            quadtree,
            camera: camera.clone(),
        }
    }
}

pub struct RenderProgress {
    pub loops: usize,
    pub raw_image: RawImage,
}

#[derive(Debug)]
pub enum RenderCommand {
    StopRender,
}

pub trait Renderer {
    fn endless_render(
        &self,
        world: World,
        tx: SyncSender<RenderProgress>,
        rx: Receiver<RenderCommand>,
    );
}
