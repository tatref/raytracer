#![allow(unused)]

use std::f64::consts::PI;

use glam::{DVec2, IVec2};
use rand::{Rng, random_range, rng, seq::IndexedRandom};
use raytracer::{
    Color,
    img::{Blending, RawImage, ToneMappingMethod},
    librt2d_forward::*,
    rt_common::*,
    worlds::*,
};

fn main() {
    let width = 800;
    let height = 600;
    let recursion_limit = 2;
    let lambda_samples = 1;

    let render_params = ForwardRenderParams {
        width,
        height,
        stop_condition: StopCondition::MaxLoops(1),
        recursion_limit,
        lambda_samples,
    };

    let max = 60;
    let chrono = std::time::Instant::now();

    let camera = Camera {
        center: DVec2::new(width as f64, height as f64) / 2.,
        size: DVec2::new(width as f64, height as f64) / 2.,
    };
    let world = simple_world(&camera, 0., 0);
    let forward_renderer = ForwardRenderer::new(render_params);

    let mut raw_image = forward_renderer.global_render(&world);
    annotate(&mut raw_image, &camera, DVec2::ZERO);

    let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard { param: 1. });
    image.save("out/out.png").unwrap();
    raw_image.to_heatmap().save("out/heatmap.png").unwrap();

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
