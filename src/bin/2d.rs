use std::{collections::HashSet, f64::consts::PI};

use glam::{DVec2, IVec2, Vec3};
use itertools::{Itertools, iproduct};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use raytracer::{
    img::{Blending, RawImage, ToneMappingMethod},
    rt::Color,
};

#[derive(Clone, Copy)]
struct Circle {
    center: DVec2,
    r: f64,
}

impl Circle {
    fn new(center: DVec2, r: f64) -> Self {
        Circle { center, r }
    }

    fn is_inside(&self, p: DVec2) -> bool {
        (self.center - p).length_squared() < self.r.powi(2)
    }

    fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
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

#[derive(Clone, Copy)]
struct Bezier {
    a: DVec2,
    b: DVec2,
    c: DVec2,
}

impl Bezier {
    fn new(a: DVec2, b: DVec2, c: DVec2) -> Self {
        Self { a, b, c }
    }

    fn at(&self, t: f64) -> DVec2 {
        let ab = self.a.lerp(self.b, t);
        let bc = self.b.lerp(self.c, t);
        ab.lerp(bc, t)
    }

    fn resolution(&self, resolution: usize) -> Vec<Segment> {
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

#[derive(Clone, Copy)]
struct Segment {
    a: DVec2,
    b: DVec2,
    n: DVec2,
}

impl Segment {
    fn new(a: DVec2, b: DVec2) -> Self {
        let n = (b - a).perp().normalize();
        Segment { a, b, n }
    }

    fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
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

#[derive(Clone)]
struct Strip {
    strip: Vec<Segment>,
}

impl Strip {
    fn hit(&self, ray: &Ray2d) -> Option<Hit2d> {
        let mut hits: Vec<Hit2d> = self
            .strip
            .iter()
            .filter_map(|segment| segment.hit(ray))
            .collect();

        hits.sort_by(|a, b| a.t.total_cmp(&b.t));

        hits.first().cloned()
    }
}

#[derive(Clone)]
enum Shape {
    Circle(Circle),
    Segment(Segment),
    Strip(Strip),
    //Bezier(Bezier),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Inside,
    Outside,
}

#[derive(Copy, Clone)]
struct Hit2d {
    t: f64,
    p: DVec2,
    n: DVec2,
    side: Side,
}

struct Ray2d {
    origin: DVec2,
    dir: DVec2,
}

impl Ray2d {
    fn new(origin: DVec2, dir: DVec2) -> Self {
        Self { origin, dir }
    }

    fn rand(origin: DVec2) -> Self {
        let angle: f64 = rand::random();
        let dir = DVec2::from_angle(2. * PI * angle);

        Self { origin, dir }
    }

    fn rand_hemisphere(origin: DVec2, n: DVec2) -> Self {
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

    // impl Ray2d
    fn hit(&self, shape: &Shape) -> Option<Hit2d> {
        match shape {
            Shape::Circle(circle) => circle.hit(self),
            Shape::Segment(segment) => segment.hit(self),
            Shape::Strip(strip) => strip.hit(self),
            //Shape::Bezier(bezier) => bezier.hit(self),
        }
    }
}

#[derive(Clone, Copy)]
enum Material {
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
    Dielectric {
        ior: f64,
    },
}

impl Material {
    fn emissive_at(d: f64, color: Color) -> Self {
        let k = d as f32 * color;
        Material::Emissive {
            emission_color: k,
            d,
        }
    }
    fn directional_emissive_at(d: f64, angle: f64, color: Color) -> Self {
        let k = d as f32 * color;
        Material::DirectionalEmissive {
            emission_color: k,
            d,
            angle,
        }
    }

    fn diffuse(color: Color) -> Self {
        Material::Diffuse {
            absorption: color.clamp(Color::ZERO, Color::ONE),
        }
    }
}

#[derive(Clone)]
struct Object {
    shape: Shape,
    mat: Material,
}

impl Object {
    fn new(shape: Shape, mat: Material) -> Self {
        Self { shape, mat }
    }

    fn segment(a: DVec2, b: DVec2, mat: Material) -> Self {
        let shape = Shape::Segment(Segment::new(a, b));
        Self { shape, mat }
    }

    fn bezier(a: DVec2, b: DVec2, c: DVec2, resolution: usize, mat: Material) -> Self {
        let strip = Bezier::new(a, b, c).resolution(resolution);
        let strip = Strip { strip };

        Object {
            shape: Shape::Strip(strip),
            mat,
        }
    }
}

struct RenderParams {
    width: i32,
    height: i32,
    spp: usize,
    recursion_limit: usize,
    denoiser: Option<Denoiser>,
}

struct Denoiser {
    top: f32,
    mask_size: i32,
    resample: f32,
}

struct World {
    //light: Light,
    objects: Vec<Object>,
}

impl World {
    fn new(objects: Vec<Object>) -> Self {
        Self { objects }
    }

    fn compute_pixel(&self, pixel: IVec2, spp: usize, recursion_limit: usize) -> Color {
        let mut color = Color::ZERO;

        for _ in 0..spp {
            let p = pixel.as_dvec2() + rand::random::<DVec2>() - DVec2::ONE / 2.;
            let ray = Ray2d::rand(p);

            let col = self.trace_ray(&ray, recursion_limit);
            color += col;
        }
        color / spp as f32
    }

    fn trace_ray(&self, ray: &Ray2d, depth: usize) -> Color {
        if depth == 0 {
            return Color::ZERO;
        }

        let mut hits: Vec<(Object, Hit2d)> = self
            .objects
            .iter()
            .filter_map(|obj| ray.hit(&obj.shape).map(|hit| (obj.clone(), hit)))
            .collect();
        hits.sort_by(|(_, a), (_, b)| a.t.total_cmp(&b.t));

        //let hit = hits.first();
        let Some((obj, hit)) = hits.first() else {
            return Color::ZERO;
        };

        let col = if hit.side == Side::Inside {
            Color::ZERO
        } else {
            match obj.mat {
                Material::Diffuse { absorption } => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::rand_hemisphere(p, hit.n);

                    let col = self.trace_ray(&r, depth - 1);

                    absorption * col * 2.
                }
                Material::Emissive { emission_color, d } => {
                    let distance_coeff = 1. / (hit.t + d);
                    emission_color * distance_coeff as f32
                    //(Color::ONE + hit.n.as_vec2().extend(0.)) * distance_coeff as f32
                }
                Material::DirectionalEmissive {
                    emission_color,
                    angle,
                    d,
                } => {
                    if ray.dir.dot(-hit.n) >= angle {
                        let distance_coeff = 1. / (hit.t + d);
                        emission_color * distance_coeff as f32
                    } else {
                        Color::ZERO
                    }
                    //(Color::ONE + hit.n.as_vec2().extend(0.)) * distance_coeff as f32
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, ray.dir.reflect(hit.n));

                    self.trace_ray(&r, depth - 1)
                }
                _ => unimplemented!(),
            }
        };

        col
    }

    fn render(&self, render_params: &RenderParams) -> RawImage {
        let chrono = std::time::Instant::now();
        let mut raw_image = RawImage::new(render_params.width, render_params.height);

        println!("Initial render...");
        // parallel v1
        for i in 0..render_params.width {
            let pixels: Vec<(i32, Color)> = (0..render_params.height)
                .into_par_iter()
                .map(|j| {
                    let pixel = IVec2::new(i, j);
                    let color =
                        self.compute_pixel(pixel, render_params.spp, render_params.recursion_limit);
                    (j, color)
                })
                .collect();
            for (j, color) in pixels {
                let pixel = IVec2::new(i, j);
                raw_image
                    .draw_pixel(pixel, color, Blending::Replace)
                    .unwrap();
            }
        }

        let elapsed = chrono.elapsed();
        dbg!(elapsed);

        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
        image.save("raw.png").unwrap();

        let Some(denoiser) = &render_params.denoiser else {
            return raw_image;
        };

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

        let chrono = std::time::Instant::now();
        println!("Removing noise...");

        let xxx: Vec<_> = recompute_pixels
            .par_iter()
            .map(|&pixel| {
                let color = self.compute_pixel(
                    pixel,
                    (render_params.spp as f32 * denoiser.resample) as usize,
                    render_params.recursion_limit,
                );
                (pixel, color)
            })
            .collect();
        for (pixel, color) in xxx {
            let _ = raw_image.draw_pixel(pixel, color, Blending::Replace);
        }

        let elapsed = chrono.elapsed();
        dbg!(elapsed);

        raw_image
    }
}

fn main() {
    let width = 800;
    let height = 600;

    let mut objects = Vec::new();

    let light = Object::new(
        Shape::Circle(Circle::new(DVec2::new(100., 100.), 50.)),
        Material::emissive_at(50., Color::ONE * 10.),
    );
    objects.push(light);

    //let absorb = Object::new(
    //    Shape::Circle(Circle::new(DVec2::new(200., 300.), 20.)),
    //    Material::diffuse(Color::new(0.9, 0.1, 0.1)),
    //);
    //objects.push(absorb);

    let mirror = Object::new(
        Shape::Circle(Circle::new(DVec2::new(700., 100.), 40.)),
        //Material::diffuse(Color::new(1., 0., 0.)),
        Material::Reflective,
    );
    objects.push(mirror);

    //let n = 6;
    //for idx in 0..n {
    //    let t = idx as f32 / n as f32;
    //    let color = Color::new(0., 1., 0.);

    //    let reflect = Object::new(
    //        Shape::Circle(Circle::new(DVec2::new(400., idx as f64 * 40.), 10.)),
    //        Material::diffuse(color),
    //    );
    //    objects.push(reflect);
    //}

    let segment = Object::segment(
        DVec2::new(400., 200.),
        DVec2::new(300., 200.),
        Material::diffuse(Color::new(0., 0., 1.)),
    );
    objects.push(segment);

    let segment = Object::segment(
        DVec2::new(300., 200.),
        DVec2::new(300., 300.),
        Material::diffuse(Color::new(1., 1., 0.)),
    );
    objects.push(segment);

    // laser
    //let segment = Segment::new(DVec2::new(100., 500.), DVec2::new(100., 490.));
    //let segment = Object::new(
    //    Shape::Segment(segment),
    //    Material::directional_emissive_at(1000., 0.99999, Color::new(0., 0.9, 0.1) * 1000.),
    //);
    //objects.push(segment);

    //let segment = Segment::new(DVec2::new(100., 500.), DVec2::new(600., 500.));
    //let segment = Object::new(Shape::Segment(segment), Material::diffuse(Color::ONE));
    //objects.push(segment);

    //let segment = Segment::new(DVec2::new(600., 490.), DVec2::new(100., 490.));
    //let segment = Object::new(Shape::Segment(segment), Material::diffuse(Color::ONE));
    //objects.push(segment);

    let segment = Segment::new(DVec2::new(720., 480.), DVec2::new(700., 500.));
    let mirror = Object::new(Shape::Segment(segment), Material::Reflective);
    objects.push(mirror);

    let bezier = Object::bezier(
        DVec2::new(200., 50.),
        DVec2::new(300., 50.),
        DVec2::new(300., 150.),
        4,
        Material::diffuse(Color::new(1., 0., 0.)),
        //Material::Reflective,
    );
    objects.push(bezier);

    let world = World::new(objects);
    let spp = 200;
    let recursion_limit = 10;

    let denoiser = Denoiser {
        top: 0.01,
        resample: 10.,
        mask_size: 2,
    };
    let render_params = RenderParams {
        height,
        spp,
        width,
        recursion_limit,
        denoiser: Some(denoiser),
    };
    let mut raw_image = world.render(&render_params);

    for i in 0..(render_params.width / 100) {
        for j in 0..5 {
            let pixel = IVec2::new(i * 100, j);
            let _ = raw_image.draw_pixel(pixel, Color::new(1., 0., 0.), Blending::Replace);
        }
    }

    for j in 0..(render_params.height / 100) {
        for i in 0..5 {
            let pixel = IVec2::new(i, j * 100);
            let _ = raw_image.draw_pixel(pixel, Color::new(1., 0., 0.), Blending::Replace);
        }
    }

    let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    image.save("out.png").unwrap();
}
