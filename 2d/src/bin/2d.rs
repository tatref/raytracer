#![allow(unused)]

use std::f64::consts::PI;

use glam::{DVec2, IVec2};
use rand::{Rng, random_range, rng, seq::IndexedRandom};
use raytracer::{
    Color,
    img::{Blending, RawImage, ToneMappingMethod},
    librt2d::*,
    worlds::*,
};

fn main() {
    let width = 800;
    let height = 600;
    let spp = 500;
    let recursion_limit = 10;

    let denoiser = Some(Denoiser {
        top: 0.01,
        oversample_factor: 10.,
        mask_size: 2,
    });
    let render_params = RenderParams {
        height,
        spp,
        width,
        recursion_limit,
        //denoiser,
        denoiser: None,
    };

    let max = 60;
    let chrono = std::time::Instant::now();
    for idx in 0..max {
        println!("{}/{}", idx, max - 1);

        let t = idx as f64 / max as f64;
        let world = sample_world(render_params.clone(), t, idx);

        let mut raw_image = world.render();
        annotate(&mut raw_image, &world, DVec2::ZERO);

        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
        image.save(&format!("out/out_{:04}.png", idx)).unwrap();
    }
    let elapsed = chrono.elapsed();
    dbg!(elapsed);

    //let t = 0.;
    //let world = sample_world(render_params.clone(), t, 0);
    //let chrono = std::time::Instant::now();
    //let mut raw_image = world.render();
    //let elapsed = chrono.elapsed();
    //dbg!(elapsed);

    //let p = DVec2::new(
    //    render_params.width as f64 / 2.,
    //    render_params.height as f64 / 2.,
    //);

    //annotate(&mut raw_image, &world, p);
    //let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    //image.save(&format!("out.png")).unwrap();
}
