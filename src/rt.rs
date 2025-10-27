use glam::{DVec2, DVec3, UVec2, Vec3};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::img::{Blending, RawImage};

pub type Color = Vec3;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ray2d {
    origin: DVec3,
    dir: DVec3,
}

impl Ray2d {
    pub fn new(origin: DVec3, dir: DVec3) -> Self {
        Self { origin, dir }
    }
    pub fn at(&self, t: f64) -> DVec3 {
        self.origin + t * self.dir
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum HitSide {
    Inside,
    Outside,
}

pub struct Hit {
    t: f64,
    p: DVec3,
    n: DVec3,
    hit_side: HitSide,
}

pub fn random_on_hemisphere(n: DVec3) -> DVec3 {
    let v: DVec3 = rand::random::<DVec3>() - 0.5 * DVec3::ONE;
    let v = v.normalize();

    if n.dot(v) > 0. { v } else { -v }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Object {
    Sphere(Sphere),
}

impl Hittable for Object {
    fn hit(&self, ray: &Ray2d) -> Option<Hit> {
        match self {
            Object::Sphere(s) => s.hit(ray),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Sphere {
    center: DVec3,
    radius: f64,
}

impl Sphere {
    pub fn new(center: DVec3, radius: f64) -> Self {
        Self { center, radius }
    }
}
impl Hittable for Sphere {
    fn hit(&self, r: &Ray2d) -> Option<Hit> {
        let oc = self.center - r.origin;
        let a = r.dir.length_squared();
        let h = r.dir.dot(oc);
        let c = oc.length_squared() - self.radius.powi(2);

        let delta = h * h - a * c;

        if delta < 0. {
            // No intersection
            return None;
        }

        // check first root
        let t = (h - delta.sqrt()) / a;

        if t >= 0. {
            let p = r.origin + r.dir * t;
            let n = (p - self.center).normalize();

            return Some(Hit {
                t,
                p,
                n,
                hit_side: HitSide::Outside,
            });
        }

        // check second root
        let t = (h + delta.sqrt()) / a;

        if t >= 0. {
            let p = r.origin + r.dir * t;
            let n = (p - self.center).normalize();

            return Some(Hit {
                t,
                p,
                n,
                hit_side: HitSide::Inside,
            });
        }

        // both roots have t<0
        return None;
    }
}

pub trait Hittable: Sync + Send {
    fn hit(&self, ray: &Ray2d) -> Option<Hit>;
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Camera {
    position: DVec3,
    dir: DVec3,
    up: DVec3,

    focal_length: f64,
    image_width: u32,
    image_height: u32,

    pixel00: DVec3,
    pixel_delta_u: DVec3,
    pixel_delta_v: DVec3,
}

impl Camera {
    pub fn new(
        position: DVec3,
        dir: DVec3,
        up: DVec3,
        focal_length: f64,
        image_width: u32,
        image_height: u32,
    ) -> Self {
        assert!(dir.is_normalized());
        assert!(up.is_normalized());

        let viewport_height = 2.;
        let viewport_width = viewport_height * (image_width as f64 / image_height as f64);

        let right = dir.cross(up);
        let viewport_u = right * viewport_width;
        let viewport_v = up * viewport_height;

        let pixel_delta_u = viewport_u / image_width as f64;
        let pixel_delta_v = viewport_v / image_height as f64;

        let viewport_upper_left =
            position + (focal_length * dir) - viewport_u / 2. - viewport_v / 2.;
        let pixel00 = viewport_upper_left + (pixel_delta_u + pixel_delta_v) / 2.;

        Self {
            dir,
            position,
            up,

            focal_length,
            image_height,
            image_width,

            pixel00,
            pixel_delta_u,
            pixel_delta_v,
        }
    }

    pub fn get_ray(&self, i: u32, j: u32) -> Ray2d {
        let rand_pixel: DVec2 = rand::random::<DVec2>() - DVec2::ONE / 2.;

        let pixel_center = self.pixel00
            + ((i as f64 + rand_pixel.x) * self.pixel_delta_u)
            + ((j as f64 + rand_pixel.y) * self.pixel_delta_v);

        let dir = (pixel_center - self.position).try_normalize().unwrap();
        let ray = Ray2d::new(self.position, dir);

        ray
    }
}

pub struct RenderParams {
    pub spp: usize,
    pub bounce_limit: usize,
}

pub struct World {
    camera: Camera,
    objects: Vec<Object>,

    render_params: RenderParams,
}

impl World {
    pub fn new(camera: Camera, objects: Vec<Object>, render_params: RenderParams) -> Self {
        Self {
            camera,
            objects: objects.iter().copied().collect(),

            render_params,
        }
    }

    pub fn trace_ray(&self, objects: &[Object], ray: Ray2d, bounce_limit: usize) -> Color {
        if bounce_limit <= 0 {
            return Color::ZERO;
        }

        let mut hits: Vec<Hit> = objects
            .iter()
            .filter_map(|object| object.hit(&ray))
            .collect();
        hits.sort_by(|a, b| a.t.total_cmp(&b.t));

        match hits.first() {
            Some(hit) => {
                let bounce_dir = (random_on_hemisphere(hit.n) + hit.n).normalize();
                let p = hit.p + hit.n * (100. * f64::EPSILON);
                let bounce_ray = Ray2d::new(p, bounce_dir);

                0.5 * self.trace_ray(objects, bounce_ray, bounce_limit - 1)
            }
            None => Self::get_background_color(&ray),
        }
    }

    pub fn compute_pixel(&self, i: u32, j: u32) -> Color {
        let color: Color = (0..self.render_params.spp)
            .into_par_iter()
            .map(|_| {
                let ray = self.camera.get_ray(i, j);
                let sub_color = self.trace_ray(&self.objects, ray, self.render_params.bounce_limit);
                sub_color
            })
            .sum();

        color
    }

    //fn render(&self, spheres: &[Box<dyn Hittable>]) -> Img {
    pub fn render(&self) -> RawImage {
        let mut img = RawImage::new(self.camera.image_width, self.camera.image_height);

        for i in 0..self.camera.image_width {
            for j in 0..self.camera.image_height {
                if i < 10 && j < 10 {
                    continue;
                }

                let color = self.compute_pixel(i, j);

                img.draw_pixel(
                    UVec2::new(i, self.camera.image_height - j - 1),
                    color,
                    &Blending::Replace,
                );
            }
        }

        img
    }

    pub fn get_background_color(r: &Ray2d) -> Color {
        let a = 0.5 * (r.dir.y as f32 + 1.);
        let color = (1. - a) * Color::ONE + a * Color::new(0., 0.7, 1.);

        let sun = true;

        let sun_dir = DVec3::new(0.3, 1., 0.).normalize();
        if sun && r.dir.dot(sun_dir) > 0.99 {
            Color::new(240., 255., 170.)
        } else {
            color
        }
    }
}
