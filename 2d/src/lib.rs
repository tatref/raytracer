#![allow(dead_code)]
#![allow(unused)]

use glam::Vec3;

pub type Color = Vec3;

pub mod aabb;
pub mod distributed;
pub mod img;
pub mod librt2d;
pub mod worlds;

pub mod spectrum;

mod noise_loop {
    use std::f64::consts::PI;

    use glam::DVec2;
    use noise_functions::{CellDistance, Noise, NoiseFn, OpenSimplex2s, OpenSimplexNoise, Perlin};
    pub struct NoiseLoop {
        center: DVec2,
        r: f64,
        zoom: f64,
        scale: f64,
    }
    impl NoiseLoop {
        pub fn new(center: DVec2, r: f64, zoom: f64, scale: f64) -> Self {
            Self {
                center,
                r,
                zoom,
                scale,
            }
        }
        /// t: [0, 1]
        pub fn at(&self, t: f64, z: f64) -> f64 {
            let angle = 2. * PI * t;
            let p = self.center + DVec2::from_angle(angle) * self.r;
            let p = p.extend(z) * self.zoom;

            let val = Perlin.add_seed(0).sample3(p.as_vec3());
            val as f64 * self.scale
        }
    }
}
