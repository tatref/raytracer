#![allow(dead_code)]
#![allow(unused_imports)]

//
// https://raytracing.github.io/books/RayTracingInOneWeekend.html
// https://bheisler.github.io/post/writing-raytracer-in-rust-part-3/
//

use std::f64::consts::PI;

use glam::{DVec2, DVec3, UVec2, Vec2, Vec3};
use image::{ImageBuffer, Rgba};

use rand::Rng;
use rayon::prelude::*;
use raytracer::{img::*, rt::*};
use serde::{Deserialize, Serialize};

fn main() {
    for cam_angle in [PI / 16.] {
        for focal_length in [2.] {
            let camera_dir = -DVec3::Z.rotate_towards(DVec3::X, cam_angle);
            let camera = Camera::new(
                DVec3::ZERO + DVec3::Y * 2.,
                camera_dir,
                DVec3::Y,
                focal_length,
                640,
                480,
            );

            let mut objects: Vec<Object> = Vec::new();
            for i in -3..3 {
                for j in -3..3 {
                    for k in -3..3 {
                        let p = i as f64 * 3. * DVec3::X + j as f64 * 3. * DVec3::Y
                            - 20. * DVec3::Z
                            + k as f64 * 3. * DVec3::Z;
                        let s = Sphere::new(p, 0.8);
                        objects.push(Object::Sphere(s));
                    }
                }
            }
            let s = Sphere::new(DVec3::new(0., 5., -20.), 5.);
            objects.push(Object::Sphere(s));

            let ground: Sphere = Sphere::new(DVec3::new(0., -10000., -20.), 10000.);
            objects.push(Object::Sphere(ground));

            let render_params = RenderParams {
                spp: 32,
                bounce_limit: 10,
            };
            let world = World::new(camera, objects, render_params);

            let chrono = std::time::Instant::now();
            println!("Rendering...");
            let raw_image = world.render();
            let elapsed = chrono.elapsed();

            println!("{:?}", elapsed);

            for (idx, tone_mapping) in [
                ToneMappingMethod::Reinhard,
                ToneMappingMethod::Gamma { gamma: 2.2 },
                ToneMappingMethod::Gamma { gamma: 2. },
                ToneMappingMethod::Gamma { gamma: 2.4 },
            ]
            .iter()
            .enumerate()
            {
                let image = raw_image.convert_to_image(tone_mapping);
                image
                    .save(&format!("out_foc={}_idx={}.png", focal_length, idx))
                    .unwrap();
            }
        }
    }
}
