use std::f64::consts::PI;

use glam::{DVec2, DVec3, UVec2};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
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
}

#[derive(Clone, Copy)]
enum Shape {
    Circle(Circle),
    //Segment(Segment)
}

#[derive(PartialEq, Eq)]
enum Side {
    Inside,
    Outside,
}

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

    fn hit(&self, shape: &Shape) -> Option<Hit2d> {
        match shape {
            Shape::Circle(circle) => {
                let a = 1.;
                let b = 2. * self.dir.dot(self.origin - circle.center);
                let c = (self.origin - circle.center).length_squared() - circle.r.powi(2);
                let delta = b.powi(2) - 4. * a * c;

                if delta < 0. {
                    return None;
                }

                let t0 = (-b - delta.sqrt()) / (2. * a);
                if t0 > 0. {
                    let p = self.origin + t0 * self.dir;
                    let n = (p - circle.center).normalize();
                    let side = if circle.is_inside(self.origin) {
                        Side::Inside
                    } else {
                        Side::Outside
                    };

                    let hit = Hit2d { t: t0, p, n, side };
                    Some(hit)
                } else {
                    let t1 = (-b + delta.sqrt()) / (2. * a);
                    if t1 > 0. {
                        let p = self.origin + t1 * self.dir;
                        let n = (p - circle.center).normalize();
                        let side = if circle.is_inside(self.origin) {
                            Side::Inside
                        } else {
                            Side::Outside
                        };

                        let hit = Hit2d { t: t1, p, n, side };
                        Some(hit)
                    } else {
                        None
                    }
                }
            }
            _ => unimplemented!(),
        }
    }
}

#[derive(Clone, Copy)]
enum Material {
    Emissive { emission_color: Color },
    Diffuse { absorption: Color },
    Reflective,
}

impl Material {
    fn emissive_at(d: f64, color: Color) -> Self {
        let k = d as f32 * color;
        Material::Emissive { emission_color: k }
    }

    fn diffuse(color: Color) -> Self {
        Material::Diffuse {
            absorption: color.clamp(Color::ZERO, Color::ONE),
        }
    }
}

#[derive(Clone, Copy)]
struct Object {
    shape: Shape,
    mat: Material,
}

impl Object {
    fn new(shape: Shape, mat: Material) -> Self {
        Self { shape, mat }
    }
}

struct RenderParams {
    width: u32,
    height: u32,
    spp: usize,
    recursion_limit: usize,
}

struct World {
    //light: Light,
    objects: Vec<Object>,
}

impl World {
    fn new(objects: Vec<Object>) -> Self {
        Self { objects }
    }

    fn compute_pixel(&self, pixel: UVec2, render_params: &RenderParams) -> Color {
        let mut color = Color::ZERO;

        for _ in 0..render_params.spp {
            let p = pixel.as_dvec2() + rand::random::<DVec2>() - DVec2::ONE / 2.;
            let ray = Ray2d::rand(p);

            let col = self.trace_ray(&ray, render_params.recursion_limit);
            color += col;
        }
        color
    }

    fn trace_ray(&self, ray: &Ray2d, depth: usize) -> Color {
        if depth == 0 {
            return Color::ZERO;
        }

        let mut hits: Vec<(Object, Hit2d)> = self
            .objects
            .iter()
            .filter_map(|&obj| ray.hit(&obj.shape).map(|hit| (obj.clone(), hit)))
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
                Material::Emissive { emission_color } => {
                    let distance_coeff = 1. / hit.t;
                    //let distance_coeff = 1. / (hit.t + obj.shape.r);
                    emission_color * distance_coeff as f32
                    //(Color::ONE + hit.n.as_vec2().extend(0.)) * distance_coeff as f32
                }
                Material::Reflective => {
                    // recurse
                    let p = hit.p + hit.n * 10000. * f64::EPSILON;
                    let r = Ray2d::new(p, hit.n);

                    self.trace_ray(&r, depth - 1)
                }
            }
        };

        col
    }

    fn render(&self, render_params: RenderParams) -> RawImage {
        let chrono = std::time::Instant::now();
        let mut raw_image = RawImage::new(render_params.width, render_params.height);

        // interative version
        //for i in 0..render_params.width {
        //    for j in 0..render_params.height {
        //        let pixel = UVec2::new(i, j);
        //        let color = self.compute_pixel(pixel, &render_params);
        //        raw_image.draw_pixel(pixel, color, &Blending::Replace);
        //    }
        //}

        // parallel v1
        for i in 0..render_params.width {
            let pixels: Vec<(u32, Color)> = (0..render_params.height)
                .into_par_iter()
                .map(|j| {
                    let pixel = UVec2::new(i, j);
                    let color = self.compute_pixel(pixel, &render_params);
                    (j, color)
                })
                .collect();
            for (j, color) in pixels {
                let pixel = UVec2::new(i, j);
                raw_image.draw_pixel(pixel, color / render_params.spp as f32, &Blending::Replace);
            }
        }

        // Parallel v2
        //let pixels: Vec<(u32, Vec<Color>)> = (0..render_params.width)
        //    .into_par_iter()
        //    .map(|i| {
        //        // test
        //        let pixels_col: Vec<Color> = (0..render_params.height)
        //            .map(|j| {
        //                let pixel = UVec2::new(i, j);
        //                self.compute_pixel(pixel, &render_params)
        //            })
        //            .collect();
        //        (i, pixels_col)
        //    })
        //    .collect();
        //for (i, pixels_col) in pixels {
        //    for (j, &color) in pixels_col.iter().enumerate() {
        //        let pixel = UVec2::new(i, j as u32);
        //        raw_image.draw_pixel(pixel, color, &Blending::Replace);
        //    }
        //}

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

    let absorb = Object::new(
        Shape::Circle(Circle::new(DVec2::new(200., 300.), 20.)),
        Material::diffuse(Color::new(0.9, 0.1, 0.1)),
    );
    objects.push(absorb);

    let reflect = Object::new(
        Shape::Circle(Circle::new(DVec2::new(500., 200.), 40.)),
        Material::Reflective,
    );
    objects.push(reflect);

    let n = 6;
    for idx in 0..n {
        let t = idx as f32 / n as f32;
        let color = Color::new(0., 1., 0.);

        let reflect = Object::new(
            Shape::Circle(Circle::new(DVec2::new(400., idx as f64 * 40.), 10.)),
            Material::diffuse(color),
        );
        objects.push(reflect);
    }

    let world = World::new(objects);
    let spp = 200;
    let recursion_limit = 5;
    let render_params = RenderParams {
        height,
        spp,
        width,
        recursion_limit,
    };
    let raw_image = world.render(render_params);

    let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    image.save("out.png").unwrap();
}
