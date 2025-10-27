#![allow(dead_code)]
#![allow(unused_imports)]

//
// https://raytracing.github.io/books/RayTracingInOneWeekend.html
//

use std::f32::consts::PI;

use glam::{UVec2, Vec3};
use image::{ImageBuffer, Rgba};

mod img;
use crate::img::*;

#[derive(Debug, Clone)]
struct Ray {
    origin: Vec3,
    dir: Vec3,
}

impl Ray {
    fn new(origin: Vec3, dir: Vec3) -> Self {
        Self { origin, dir }
    }
    fn at(&self, t: f32) -> Vec3 {
        self.origin + t * self.dir
    }
}

struct Hit {
    t: f32,
    p: Vec3,
    n: Vec3,
}

struct Sphere {
    center: Vec3,
    radius: f32,
}

impl Sphere {
    fn new(center: Vec3, radius: f32) -> Self {
        Self { center, radius }
    }
}
impl Hittable for Sphere {
    fn hit(&self, r: &Ray) -> Option<Hit> {
        let oc = self.center - r.origin;
        let a = r.dir.length_squared();
        let h = r.dir.dot(oc);
        let c = oc.length_squared() - self.radius.powi(2);

        let delta = h * h - a * c;

        if delta < 0. {
            // No intersection
            return None;
        }

        let t = (h - delta.sqrt()) / a;

        let p = r.origin + r.dir * t;
        let n = (p - self.center).normalize();

        return Some(Hit {
            t,
            p,
            n,
            hit_side: HitSide::Outside,
        });
    }
}

trait Hittable {
    fn hit(&self, ray: &Ray) -> Option<Hit>;
}

struct Camera {
    position: Vec3,
    dir: Vec3,
    up: Vec3,

    focal_length: f32,
    image_width: u32,
    image_height: u32,

    pixel00: Vec3,
    pixel_delta_u: Vec3,
    pixel_delta_v: Vec3,
}

impl Camera {
    fn new(
        position: Vec3,
        dir: Vec3,
        up: Vec3,
        focal_length: f32,
        image_width: u32,
        image_height: u32,
    ) -> Self {
        assert!(dir.is_normalized());
        assert!(up.is_normalized());

        let viewport_height = 2.;
        let viewport_width = viewport_height * (image_width as f32 / image_height as f32);

        let right = dir.cross(up);
        let viewport_u = right * viewport_width;
        let viewport_v = up * viewport_height;

        let pixel_delta_u = viewport_u / image_width as f32;
        let pixel_delta_v = viewport_v / image_height as f32;

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

    fn render(&self, spheres: &[Box<dyn Hittable>]) -> Img {
        let mut img = Img::new(self.image_width, self.image_height);

        for i in 0..self.image_width {
            for j in 0..self.image_height {
                if i < 10 && j < 10 {
                    continue;
                }

                let mut color = Vec3::ZERO;

                let spp = 2;

                let sub_pixel_delta_u = self.pixel_delta_u / spp as f32;
                let sub_pixel_delta_v = self.pixel_delta_v / spp as f32;

                let pixel_center = self.pixel00
                    + (i as f32 * self.pixel_delta_u)
                    + (j as f32 * self.pixel_delta_v);

                let dir = (pixel_center - self.position).try_normalize().unwrap();

                let ray = Ray::new(self.position, dir);

                let mut hits: Vec<Hit> = spheres
                    .iter()
                    .filter_map(|sphere| sphere.hit(&ray))
                    .collect();
                hits.sort_by(|a, b| a.t.total_cmp(&b.t));

                let sub_color = if let Some(hit) = hits.first() {
                    hit.n
                } else {
                    Self::get_background_color(&ray)
                };

                color += sub_color;
                img.draw_pixel(UVec2::new(i, self.image_height - j - 1), color);
            }
        }

        img
    }

    fn get_background_color(r: &Ray) -> Vec3 {
        let a = 0.5 * (r.dir.y + 1.);
        let color = (1. - a) * Vec3::ONE + a * Vec3::new(0., 0.7, 1.);

        color
    }
}

fn main() {
    for focal_length in [0.5, 1., 2.] {
        let camera_dir = -Vec3::Z.rotate_towards(Vec3::X, cam_angle);
        let camera = Camera::new(Vec3::ZERO, camera_dir, Vec3::Y, focal_length, 640, 480);

        let mut spheres: Vec<Box<dyn Hittable>> = Vec::new();
        for i in -3..3 {
            for j in -3..3 {
                for k in -3..3 {
                    let p = i as f32 * 3. * Vec3::X + j as f32 * 3. * Vec3::Y - 20. * Vec3::Z
                        + k as f32 * 3. * Vec3::Z;
                    let s = Sphere::new(p, 0.8);
                    spheres.push(Box::new(s));
                }
            }
        }
        //let s = Sphere::new(-Vec3::Z * 20., 5.);
        //spheres.push(Box::new(s));

        let ground: Sphere = Sphere::new(Vec3::Y * -1000. + Vec3::Z * -50., 1000.);
        spheres.push(Box::new(ground));

        let chrono = std::time::Instant::now();
        let img = camera.render(&spheres);
        let elapsed = chrono.elapsed();

        println!("{:?}", elapsed);

        let image = img.into_image();
        image
            .save(&format!("out_foc={}.png", focal_length))
            .unwrap();
    }
}
