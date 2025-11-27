#![allow(unused)]

use std::f64::consts::PI;

use glam::{DVec2, IVec2};
use rand::{Rng, random_range, rng, seq::IndexedRandom};
use raytracer::{
    Color,
    img::{Blending, RawImage, ToneMappingMethod},
    librt2d::*,
};

use noise_functions::{CellDistance, Noise, NoiseFn, OpenSimplex2s, OpenSimplexNoise, Perlin};
struct NoiseLoop {
    center: DVec2,
    r: f64,
    zoom: f64,
    scale: f64,
}
impl NoiseLoop {
    fn new(center: DVec2, r: f64, zoom: f64, scale: f64) -> Self {
        Self {
            center,
            r,
            zoom,
            scale,
        }
    }
    /// t: [0, 1]
    fn at(&self, t: f64, z: f64) -> f64 {
        let angle = 2. * PI * t;
        let p = self.center + DVec2::from_angle(angle) * self.r;
        let p = p.extend(z) * self.zoom;

        let val = Perlin.add_seed(0).sample3(p.as_vec3());
        val as f64 * self.scale
    }
}

fn sample_world(render_params: RenderParams, t: f64, idx: u64) -> World {
    let mut objects = Vec::new();

    use rand::RngCore;
    use rand::prelude::*;
    use rand_chacha::ChaCha20Rng;

    let mut rng = ChaCha20Rng::seed_from_u64(2);

    //for _ in 0..50 {
    //    let center = rand::random::<DVec2>() * DVec2::new(300., 200.) + DVec2::new(200., 300.);
    //    let light = Object::new(
    //        Shape::Circle(Circle::new(center, 1.)),
    //        Material::emissive_at(50., Color::ONE * 10.),
    //    );
    //    objects.push(light);
    //}

    let light = Object::new(
        Shape::Circle(Circle::new(DVec2::new(600., 400.), 20.)),
        Material::emissive_at(50., Color::ONE * 10.),
    );
    objects.push(light);

    let angle = 2. * PI * t * 2.;
    let phi = PI * 2. / 3.;
    let dielectric = Object::new(
        Shape::Circle(Circle::new(
            DVec2::new(
                600. + 50. * (angle + 0. * phi).sin(),
                400. + 50. * (angle + 0. * phi).cos(),
            ),
            10.,
        )),
        Material::dieletric(0.7),
    );
    objects.push(dielectric);

    let diffuse = Object::new(
        Shape::Circle(Circle::new(
            DVec2::new(
                600. + 50. * (angle + 1. * phi).sin(),
                400. + 50. * (angle + 1. * phi).cos(),
            ),
            10.,
        )),
        Material::diffuse(Color::ONE),
    );
    objects.push(diffuse);

    let phi_2 = phi / 10.;
    let a = DVec2::new(
        600. + 50. * (angle + 2. * phi - phi_2).sin(),
        400. + 50. * (angle + 2. * phi - phi_2).cos(),
    );
    let b = DVec2::new(
        600. + 50. * (angle + 2. * phi + phi_2).sin(),
        400. + 50. * (angle + 2. * phi + phi_2).cos(),
    );
    let diffuse = Object::new(Shape::Segment(Segment::new(a, b)), Material::Reflective);
    objects.push(diffuse);

    let big_light = Object::new(
        Shape::Circle(Circle::new(DVec2::new(50. + angle.sin() * 100., 100.), 20.)),
        Material::emissive_at(50., Color::ONE * 10.),
    );
    objects.push(big_light);

    // wavy bezier
    let start_x = 400.;
    let end_x = 800.;
    let y = 100.;

    let n_points = 30;
    let mut points = Vec::new();
    //points.push(DVec2::new(400., 100.));

    let noise_loop = NoiseLoop::new(DVec2::ZERO, 1., 1., 10.);

    for i in 0..n_points {
        let x = start_x + (end_x - start_x) * (i as f64 / n_points as f64);
        let dy = /*noise_loop.at(t, i as f64 / n_points as f64) */
             ((i as f64 / n_points as f64 + t) * 2. * PI).sin() * 10.;
        let p = DVec2::new(x, y + dy);
        points.push(p);
    }
    //points.push(DVec2::new(800., 100.));

    let bezier = Bezier::new(&points);
    let strip = bezier.as_segments(100);
    let bezier_strip = Strip::new(&strip);
    //let bezier_obj = Object::new(Shape::Strip(bezier_strip), Material::diffuse(Color::ONE));
    let bezier_obj = Object::new(Shape::Strip(bezier_strip), Material::dieletric(0.7));
    objects.push(bezier_obj);

    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(400., 0.), DVec2::new(800., 0.))),
        Material::diffuse(Color::ONE),
    );

    // ?
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(0., 250.), DVec2::new(100., 250.))),
        Material::diffuse(Color::ONE),
    );
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(110., 250.), DVec2::new(200., 250.))),
        Material::diffuse(Color::ONE),
    );
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(200., 250.), DVec2::new(200., 500.))),
        Material::diffuse(Color::ONE),
    );
    objects.push(wall);

    // Box
    let center = DVec2::new(350., 250.);
    let r = 50.;
    let angle = 2. * PI * t;
    let d_phi = 2. * PI / 4.;
    let segment = Object::segment(
        center + DVec2::from_angle(angle + 1. * d_phi) * r,
        center + DVec2::from_angle(angle + 0. * d_phi) * r,
        Material::diffuse(Color::new(0., 0., 1.)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 2. * d_phi) * r,
        center + DVec2::from_angle(angle + 1. * d_phi) * r,
        Material::diffuse(Color::new(1., 1., 0.)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 3. * d_phi) * r,
        center + DVec2::from_angle(angle + 2. * d_phi) * r,
        Material::diffuse(Color::new(0., 1., 0.)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 4. * d_phi) * r,
        center + DVec2::from_angle(angle + 3. * d_phi) * r,
        Material::diffuse(Color::new(0., 1., 1.)),
    );
    objects.push(segment);

    let dieletric = Object::from_points(
        &[
            DVec2::new(50., 300.),
            DVec2::new(200., 300.),
            DVec2::new(125., 350.),
        ],
        Material::dieletric(0.7),
        //Material::diffuse(Color::ZERO),
    );
    objects.push(dieletric);

    // laser
    //let segment = Segment::new(DVec2::new(220., 500.), DVec2::new(220., 490.));
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
        &[
            DVec2::new(200., 50.),
            DVec2::new(300., 50.),
            DVec2::new(300., 150.),
        ],
        4,
        //Material::diffuse(Color::new(1., 0., 0.)),
        Material::Reflective,
    );
    objects.push(bezier);

    let world = World::new(objects, render_params);

    world
}

fn simple_world(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let center = DVec2::new(400., 300.);
    let light = Object::new(
        Shape::Circle(Circle::new(center, 1.)),
        Material::emissive_at(50., Color::ONE * 20.),
    );
    objects.push(light);

    let segment = Segment::new(DVec2::new(300., 200.), DVec2::new(500., 200.));
    let segment = Object::new(Shape::Segment(segment), Material::diffuse(Color::ONE));
    objects.push(segment);

    let world = World::new(objects, render_params);

    world
}

fn complex_world(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let materials = [
        Material::Reflective,
        Material::emissive_at(10., Color::ONE * 3.),
        Material::diffuse(Color::ONE),
        Material::dieletric(0.7),
    ];
    let mut rng = rand::rng();
    for _ in 0..4000 {
        //let mat = materials.choose(&mut rng).unwrap();
        let mat = Material::emissive_at(10., Color::ONE * 3.);
        let center = rand::random::<DVec2>() * DVec2::new(800., 600.);
        let light = Object::new(Shape::Circle(Circle::new(center, 1.)), mat.clone());
        objects.push(light);
    }

    let world = World::new(objects, render_params);

    world
}

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
