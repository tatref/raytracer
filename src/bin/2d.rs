#![allow(unused)]

use std::f64::consts::PI;

use glam::{DVec2, IVec2};
use rand::seq::IndexedRandom;
use raytracer::{
    img::{Blending, RawImage, ToneMappingMethod},
    librt2d::*,
    rt::Color,
};

fn sample_world(t: f64) -> World {
    let mut objects = Vec::new();

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
        Material::dieletric(1.5),
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
        Material::dieletric(1.5),
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
        DVec2::new(200., 50.),
        DVec2::new(300., 50.),
        DVec2::new(300., 150.),
        4,
        Material::diffuse(Color::new(1., 0., 0.)),
        //Material::Reflective,
    );
    objects.push(bezier);

    let world = World::new(objects);

    world
}

fn simple_world(_t: f64) -> World {
    let mut objects = Vec::new();

    let center = rand::random::<DVec2>() * DVec2::new(300., 200.) + DVec2::new(200., 300.);
    let light = Object::new(
        Shape::Circle(Circle::new(center, 1.)),
        Material::emissive_at(50., Color::ONE * 10.),
    );
    objects.push(light);

    let world = World::new(objects);

    world
}

fn complex_world(_t: f64) -> World {
    let mut objects = Vec::new();

    let materials = [
        Material::Reflective,
        Material::emissive_at(10., Color::ONE * 3.),
        Material::diffuse(Color::ONE),
        Material::dieletric(1.5),
    ];
    let mut rng = rand::rng();
    for _ in 0..4000 {
        //let mat = materials.choose(&mut rng).unwrap();
        let mat = Material::emissive_at(10., Color::ONE * 3.);
        let center = rand::random::<DVec2>() * DVec2::new(800., 600.);
        let light = Object::new(Shape::Circle(Circle::new(center, 1.)), mat.clone());
        objects.push(light);
    }

    let world = World::new(objects);

    world
}

fn annotate(raw_image: &mut RawImage, render_params: &RenderParams) {
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
}

fn main() {
    let width = 800;
    let height = 600;
    let spp = 5000;
    let recursion_limit = 10;

    let denoiser = Denoiser {
        top: 0.001,
        oversample_factor: 10.,
        mask_size: 2,
    };
    let render_params = RenderParams {
        height,
        spp,
        width,
        recursion_limit,
        denoiser: None,
    };

    let max = 100;
    let chrono = std::time::Instant::now();
    for idx in 0..max {
        println!("{}/{}", idx, max - 1);

        let t = idx as f64 / max as f64;
        let world = sample_world(t);

        let mut raw_image = world.render(&render_params);
        annotate(&mut raw_image, &render_params);

        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
        image.save(&format!("out/out_{:04}.png", idx)).unwrap();
    }
    let elapsed = chrono.elapsed();
    dbg!(elapsed);

    let t = 0.;
    let world = sample_world(t);
    let mut raw_image = world.render(&render_params);
    annotate(&mut raw_image, &render_params);
    let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
    image.save(&format!("out.png")).unwrap();
}
