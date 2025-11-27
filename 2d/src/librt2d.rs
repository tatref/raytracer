use std::{collections::HashSet, f64::consts::PI};

use crate::{
    Color,
    img::{Blending, RawImage, ToneMappingMethod},
};
use glam::{DVec2, IVec2, Vec3};
use itertools::{Itertools, iproduct};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

pub fn annotate(raw_image: &mut RawImage, world: &World, p: DVec2) {
    for i in 0..(world.render_params.width / 100) {
        for j in 0..5 {
            let pixel = IVec2::new(i * 100, j);
            let _ = raw_image.draw_pixel(pixel, Color::new(1., 0., 0.), Blending::Replace);
        }
    }

    for j in 0..(world.render_params.height / 100) {
        for i in 0..5 {
            let pixel = IVec2::new(i, j * 100);
            let _ = raw_image.draw_pixel(pixel, Color::new(1., 0., 0.), Blending::Replace);
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

#[derive(Serialize, Deserialize)]
pub struct Node {
    objects: Vec<Object>,
    aabb: Aabb,
    top_right: Option<Box<Node>>,
    top_left: Option<Box<Node>>,
    bottom_left: Option<Box<Node>>,
    bottom_right: Option<Box<Node>>,
}

#[derive(Copy, Clone)]
pub enum Direction {
    TopRight,
    TopLeft,
    BottomLeft,
    BottomRight,
}

impl Node {
    pub fn new(aabb: Aabb) -> Self {
        let objects = Vec::new();
        Self {
            objects,
            aabb,
            top_right: None,
            top_left: None,
            bottom_left: None,
            bottom_right: None,
        }
    }

    pub fn insert(&mut self, obj: Object) -> bool {
        use Direction::*;
        let obj_aabb = obj.aabb();

        // too big to fit in children
        if obj_aabb.half_size.cmpgt(self.aabb.half_size / 2.).any() {
            self.objects.push(obj);
            return true;
        }

        let sub_aabbs = [BottomLeft, BottomRight, TopLeft, TopRight].map(|dir| self.sub_aabb(dir));

        for (node, sub_aabb) in [
            &mut self.bottom_left,
            &mut self.bottom_right,
            &mut self.top_left,
            &mut self.top_right,
        ]
        .iter_mut()
        .zip(sub_aabbs)
        {
            if sub_aabb.contains(&obj.aabb()) {
                let child = node.get_or_insert_with(|| Box::new(Node::new(sub_aabb)));
                return child.insert(obj);
            }
        }

        return false;
    }

    pub fn sub_aabb(&self, dir: Direction) -> Aabb {
        let half_size = self.aabb.half_size / 2.;
        let center = match dir {
            Direction::TopRight => self.aabb.center + DVec2::new(half_size.x, half_size.y),
            Direction::TopLeft => self.aabb.center + DVec2::new(-half_size.x, half_size.y),
            Direction::BottomLeft => self.aabb.center + DVec2::new(-half_size.x, -half_size.y),
            Direction::BottomRight => self.aabb.center + DVec2::new(half_size.x, -half_size.y),
        };
        Aabb::new(center, half_size)
    }

    pub fn hit(&self, ray: &Ray2d) -> Option<(Object, Hit2d)> {
        if self.aabb.hit(&ray).is_none() {
            // ray does not hit current node nor children
            return None;
        }

        // check if ray hits children
        let children_hits: Vec<(Object, Hit2d)> = [
            &self.bottom_left,
            &self.bottom_right,
            &self.top_left,
            &self.top_right,
        ]
        .iter()
        .filter_map(|node| {
            node.as_ref().map(|node| {
                //
                match node.aabb.hit(ray) {
                    None => None,
                    Some(_) => node.hit(ray), // recurse
                }
            })
        })
        .flatten()
        .collect();

        // object could be in current node?
        let local_hits: Vec<(Object, Hit2d)> = self
            .objects
            .iter()
            .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj.clone(), hit)))
            .collect();

        let mut hits = Vec::new();
        hits.extend(children_hits);
        hits.extend(local_hits);

        hits.sort_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));

        hits.first().cloned()
    }

    pub fn draw(&self, raw_image: &mut RawImage, color: Color) {
        self.aabb.draw(raw_image, color);

        for sub_node in [
            &self.top_right,
            &self.top_left,
            &self.bottom_right,
            &self.bottom_left,
        ] {
            match sub_node {
                Some(sub_node) => sub_node.draw(raw_image, color),
                None => (),
            }
        }
    }
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct Aabb {
    center: DVec2,
    half_size: DVec2,
}
impl Aabb {
    pub fn new(center: DVec2, half_size: DVec2) -> Self {
        Self { center, half_size }
    }

    pub fn contains(&self, other: &Self) -> bool {
        (self.center - self.half_size)
            .cmplt(other.center - other.half_size)
            .all()
            && (self.center + self.half_size)
                .cmpgt(other.center + other.half_size)
                .all()
    }

    pub fn union(&self, other: &Self) -> Self {
        let min = (self.center - self.half_size).min(other.center - other.half_size);
        let max = (self.center - self.half_size).max(other.center - other.half_size);

        let center = min.midpoint(max);
        let half_size = (max - min) / 2.;

        Aabb::new(center, half_size)
    }

    pub fn hit(&self, r: &Ray2d) -> Option<f64> {
        // https://tavianator.com/2011/ray_box.html
        let tx1 = (self.center.x - self.half_size.x - r.origin.x) * r.dir.recip().x;
        let tx2 = (self.center.x + self.half_size.x - r.origin.x) * r.dir.recip().x;

        let tmin = tx1.min(tx2);
        let tmax = tx1.max(tx2);

        let ty1 = (self.center.y - self.half_size.y - r.origin.y) * r.dir.recip().y;
        let ty2 = (self.center.y + self.half_size.y - r.origin.y) * r.dir.recip().y;

        let tmin = ty1.min(ty2).max(tmin);
        let tmax = ty1.max(ty2).min(tmax);

        if tmax > tmin && tmax > 0. {
            Some(tmax)
        } else {
            None
        }
    }

    pub fn draw(&self, raw_image: &mut RawImage, color: Color) {
        let blend = Blending::Replace;

        //// cross
        //for i in -3..=3 {
        //    let pixel = self.center.as_ivec2() + IVec2::new(i, 0);
        //    let _ = raw_image.draw_pixel(pixel, color, blend);
        //}
        //for j in -3..=3 {
        //    let pixel = self.center.as_ivec2() + IVec2::new(0, j);
        //    let _ = raw_image.draw_pixel(pixel, color, blend);
        //}

        // box TOP
        let h = self.half_size.x as i32;
        for i in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(i, -self.half_size.y as i32);
            let _ = raw_image.draw_pixel(pixel, color, blend);
        }
        // box BOTTOM
        for i in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(i, self.half_size.y as i32);
            let _ = raw_image.draw_pixel(pixel, color, blend);
        }
        // box LEFT
        let h = self.half_size.y as i32;
        for j in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(self.half_size.x as i32, j);
            let _ = raw_image.draw_pixel(pixel, color, blend);
        }
        // box RIGHT
        for j in -h..h {
            let pixel = self.center.as_ivec2() + IVec2::new(-self.half_size.x as i32, j);
            let _ = raw_image.draw_pixel(pixel, color, blend);
        }
    }
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

#[derive(Copy, Clone)]
pub struct Hit2d {
    t: f64,
    p: DVec2,
    n: DVec2,
    side: Side,
}

pub struct Ray2d {
    origin: DVec2,
    dir: DVec2,
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
        emission_color: Color,
        /// inner radius
        /// for better fall off
        d: f64,
    },
    DirectionalEmissive {
        /// emission color
        emission_color: Color,
        angle: f64,
        /// inner radius
        /// for better fall off
        d: f64,
    },
    Diffuse {
        absorption: Color,
    },
    Reflective,
    // FIXME: 1/ior????
    Dielectric {
        ior: f64,
    },
}

impl Material {
    pub fn emissive_at(d: f64, color: Color) -> Self {
        let k = d as f32 * color;
        Material::Emissive {
            emission_color: k,
            d,
        }
    }
    pub fn directional_emissive_at(d: f64, angle: f64, color: Color) -> Self {
        let k = d as f32 * color;
        Material::DirectionalEmissive {
            emission_color: k,
            d,
            angle,
        }
    }

    pub fn diffuse(color: Color) -> Self {
        Material::Diffuse {
            absorption: color.clamp(Color::ZERO, Color::ONE),
        }
    }

    pub fn dieletric(ior: f64) -> Self {
        Material::Dielectric { ior }
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
    pub denoiser: Option<Denoiser>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Denoiser {
    pub top: f32,
    pub mask_size: i32,
    pub oversample_factor: f32,
}

#[derive(Serialize, Deserialize)]
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

        for obj in objects.iter().cloned() {
            quadtree.insert(obj);
        }

        Self {
            objects,
            quadtree,
            render_params,
        }
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

    pub fn compute_pixel(&self, pixel: IVec2, spp: usize, recursion_limit: usize) -> Color {
        let mut color = Color::ZERO;

        for _ in 0..spp {
            let p = pixel.as_dvec2() + rand::random::<DVec2>() - DVec2::ONE / 2.;
            let ray = Ray2d::rand(p);

            let (_t, col) = self.trace_ray(&ray, recursion_limit);
            color += col;
        }
        color / spp as f32
    }

    pub fn trace_ray(&self, ray: &Ray2d, depth: usize) -> (f64, Color) {
        if depth == 0 {
            return (0., Color::ZERO);
        }

        // Vec version
        let mut hits: Vec<(Object, Hit2d)> = self
            .objects
            .iter()
            .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj.clone(), hit)))
            .collect();
        hits.sort_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));
        let Some((obj, hit)) = hits.first() else {
            return (0., Color::ZERO);
        };

        // Quadtree version
        //let hit: Option<(Object, Hit2d)> = self.quadtree.hit(ray);
        //let Some((obj, hit)) = hit else {
        //    return Color::ZERO;
        //};

        if hit.side == Side::Inside {
            match obj.mat {
                Material::Emissive { emission_color, d } => (d, emission_color),
                Material::Dielectric { ior } => {
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(-hit.n, 1. / ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(-hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let (d2, col) = self.trace_ray(&r, depth - 1);

                    (hit.t + d2, col)
                }
                _ => (0., Color::ZERO),
            }
        } else {
            // outside
            match obj.mat {
                Material::Diffuse { absorption } => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::rand_hemisphere(p, hit.n);

                    let (d2, col) = self.trace_ray(&r, depth - 1);

                    (hit.t + d2, absorption * col)
                }
                Material::Emissive { emission_color, d } => {
                    let distance_coeff = 1. / (hit.t + d);
                    (hit.t + d, emission_color * distance_coeff as f32)
                }
                Material::DirectionalEmissive {
                    emission_color,
                    angle,
                    d,
                } => {
                    if ray.dir.dot(-hit.n) >= angle {
                        let distance_coeff = 1. / (hit.t + d);
                        (hit.t + d, emission_color * distance_coeff as f32)
                    } else {
                        (0., Color::ZERO)
                    }
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, ray.dir.reflect(hit.n));

                    let (t, col) = self.trace_ray(&r, depth - 1);
                    (hit.t + t, col)
                }
                Material::Dielectric { ior } => {
                    let p = hit.p - hit.n * 10000. * f64::EPSILON;
                    let refracted_ray = ray.dir.refract(hit.n, ior);
                    let r = if refracted_ray == DVec2::ZERO {
                        Ray2d::new(p, ray.dir.reflect(hit.n))
                    } else {
                        Ray2d::new(p, refracted_ray)
                    };
                    let (t, col) = self.trace_ray(&r, depth - 1);
                    (hit.t + t, col)
                }
            }
        }
    }

    pub fn render_column(&self, render_params: &RenderParams, i: i32) -> Vec<(i32, Color)> {
        let pixels: Vec<(i32, Color)> = (0..render_params.height)
            .into_par_iter()
            .map(|j| {
                let pixel = IVec2::new(i, j);
                let color =
                    self.compute_pixel(pixel, render_params.spp, render_params.recursion_limit);
                (j, color)
            })
            .collect();
        pixels
    }

    pub fn render(&self) -> RawImage {
        let mut raw_image = RawImage::new(self.render_params.width, self.render_params.height);

        // parallel v1
        for i in 0..self.render_params.width {
            let pixels = self.render_column(&self.render_params, i);

            for (j, color) in pixels {
                let pixel = IVec2::new(i, j);
                raw_image
                    .draw_pixel(pixel, color, Blending::Replace)
                    .unwrap();
            }
        }

        //let aabb_color = Color::new(1., 0., 1.);
        //self.quadtree.draw(&mut raw_image, aabb_color);

        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
        image.save("raw.png").unwrap();

        match self.render_params.denoiser {
            Some(_) => self.denoise(&self.render_params, raw_image),
            None => raw_image,
        }
    }

    fn denoise(&self, render_params: &RenderParams, mut raw_image: RawImage) -> RawImage {
        let denoiser = render_params.denoiser.as_ref().unwrap();

        let laplace_kernel = [[-1., -1., -1.], [-1., 8., -1.], [-1., -1., -1.]];
        let laplace_image = raw_image
            .convolution(laplace_kernel)
            .map_pixel(|p| Vec3::splat(p.length()));
        laplace_image
            .convert_to_image(&ToneMappingMethod::Reinhard)
            .save("laplace.png")
            .unwrap();

        let mut sorted_pixels: Vec<(usize, &f32)> = laplace_image.data.iter().enumerate().collect();
        sorted_pixels.sort_by(|a, b| a.1.total_cmp(&b.1));
        sorted_pixels.reverse();

        let total_pixels = render_params.width * render_params.height * 3;
        let top = (total_pixels as f32 * denoiser.top) as usize;
        let recompute_indexes: Vec<usize> = sorted_pixels[0..top]
            .iter()
            .map(|(idx, _val)| *idx)
            .collect();

        let mut recompute_mask = raw_image.clone();
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
            let _ = recompute_mask.draw_pixel(*pixel, Color::new(1., 0., 1.), Blending::Replace);
        }
        recompute_mask
            .convert_to_image(&ToneMappingMethod::Reinhard)
            .save("mask.png")
            .unwrap();

        let xxx: Vec<_> = recompute_pixels
            .par_iter()
            .map(|&pixel| {
                let color = self.compute_pixel(
                    pixel,
                    (render_params.spp as f32 * denoiser.oversample_factor) as usize,
                    render_params.recursion_limit,
                );
                (pixel, color)
            })
            .collect();
        for (pixel, color) in xxx {
            let _ = raw_image.draw_pixel(pixel, color, Blending::Replace);
        }

        raw_image
    }
}
