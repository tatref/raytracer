use std::f64::consts::PI;

use glam::DVec2;
use rand::seq::IndexedRandom;

use crate::{
    Color,
    librt2d::*,
    spectrum::{Spectrum, SpectrumColor},
};

pub fn cornell_box(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let center = DVec2::new(
        render_params.width as f64 / 2.,
        render_params.height as f64 / 2.,
    );
    let light = Object::new(
        Shape::Circle(Circle::new(center, 5.)),
        Material::emissive(Spectrum::color(SpectrumColor::White) * 0.5),
    );
    objects.push(light);

    let top = Object::segment(
        DVec2::new(200., 100.),
        DVec2::new(600., 100.),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(top);

    let right = Object::segment(
        DVec2::new(600., 100.),
        DVec2::new(600., 500.),
        Material::diffuse(Spectrum::color(SpectrumColor::Red)),
    );
    objects.push(right);

    let left = Object::segment(
        DVec2::new(200., 500.),
        DVec2::new(200., 100.),
        Material::diffuse(Spectrum::color(SpectrumColor::Blue)),
    );
    objects.push(left);

    let sphere = Object::new(
        Shape::Circle(Circle::new(DVec2::new(300., 200.), 50.)),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(sphere);

    let sphere = Object::new(
        Shape::Circle(Circle::new(DVec2::new(500., 400.), 50.)),
        Material::dieletric(1.5),
    );
    objects.push(sphere);

    let sphere = Object::new(
        Shape::Circle(Circle::new(DVec2::new(300., 400.), 50.)),
        Material::Dielectric {
            //ior: Ior::Cauchy { a: 1.45, b: 0.1 },
            ior: Ior::Cauchy { a: 1.45, b: 0.05 },
        },
    );
    objects.push(sphere);

    //objects.push(sphere);
    //let sides = 3;
    //let points: Vec<DVec2> = (0..sides)
    //    .map(|i| {
    //        DVec2::from_angle(i as f64 / sides as f64 * 2. * PI + 0.5) * 50.
    //            + DVec2::new(500., 400.)
    //    })
    //    .rev()
    //    .collect();
    //let prism = Object::from_points(
    //    &points,
    //    Material::Dielectric {
    //        ior: Ior::Cauchy { a: 1.45, b: 0.04 },
    //    },
    //);
    //objects.push(prism);

    let sides = 4;
    let points: Vec<DVec2> = (0..sides)
        .map(|i| {
            DVec2::from_angle(i as f64 / sides as f64 * 2. * PI) * 50. + DVec2::new(500., 200.)
        })
        .rev()
        .collect();
    let xxx = Object::from_points(
        &points,
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(xxx);

    let world = World::new(objects, render_params);

    world
}

pub fn simple_world(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let center = DVec2::new(400., 300.);
    let light = Object::new(
        Shape::Circle(Circle::new(center, 1.)),
        Material::emissive(Spectrum::color(SpectrumColor::White) * 10.),
    );
    objects.push(light);

    let segment = Segment::new(DVec2::new(300., 200.), DVec2::new(500., 200.));
    let segment = Object::new(
        Shape::Segment(segment),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(segment);

    let world = World::new(objects, render_params);

    world
}

pub fn spectrum_world(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let spacing = 80.;

    let wall = Object::segment(
        DVec2::new(0., 300.),
        DVec2::new(800., 300.),
        Material::diffuse(Spectrum::default()),
    );
    objects.push(wall);

    for (idx, color) in SpectrumColor::iter_colors().iter().enumerate() {
        let spectrum = Spectrum::color(*color);

        // top lights
        let center = DVec2::new(100. + idx as f64 * spacing, 150.);
        let light = Object::new(
            Shape::Circle(Circle::new(center, 10.)),
            Material::emissive(spectrum * 0.1),
        );
        objects.push(light);

        // bottom lights
        let center = DVec2::new(100. + idx as f64 * spacing, 450.);
        let light = Object::new(
            Shape::Circle(Circle::new(center, 10.)),
            Material::emissive(Spectrum::color(SpectrumColor::White) * 0.1),
        );
        objects.push(light);

        // top walls
        let a = DVec2::new(100. - spacing / 2. + idx as f64 * spacing, 400.);
        let b = DVec2::new(100. - spacing / 2. + (idx + 1) as f64 * spacing, 400.);
        let mat = Material::diffuse(Spectrum::color(*color));
        //let mat = Material::diffuse(Spectrum::default());
        let colored_wall = Object::segment(a, b, mat);
        objects.push(colored_wall);

        // vertical walls
        let b = DVec2::new(100. - spacing / 2. + idx as f64 * spacing, 400.);
        let a = DVec2::new(100. - spacing / 2. + idx as f64 * spacing, 600.);
        //let mat = Material::diffuse(Spectrum::default());
        let colored_wall = Object::segment(a, b, mat);
        objects.push(colored_wall);
    }

    let world = World::new(objects, render_params);

    world
}

pub fn complex_world(render_params: RenderParams, _t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let materials = [
        Material::Reflective,
        Material::emissive(Spectrum::color(SpectrumColor::White) * 0.08),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
        Material::dieletric(1.5),
    ];
    let mut rng = rand::rng();
    for _ in 0..200 {
        let mat = materials.choose(&mut rng).unwrap();
        let center = rand::random::<DVec2>() * DVec2::new(800., 600.);
        let light = Object::new(Shape::Circle(Circle::new(center, 5.)), mat.clone());
        objects.push(light);
    }

    let world = World::new(objects, render_params);

    world
}

pub fn sample_world(render_params: RenderParams, t: f64, _idx: u64) -> World {
    let mut objects = Vec::new();

    let light = Object::new(
        Shape::Circle(Circle::new(DVec2::new(600., 400.), 20.)),
        Material::emissive(Spectrum::color(SpectrumColor::White) * 0.1),
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
        //Material::dieletric(1.5),
        Material::Dielectric {
            ior: Ior::Cauchy { a: 1.45, b: 0.0354 },
        },
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
        Material::diffuse(Spectrum::color(SpectrumColor::Red)),
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
        Material::emissive(Spectrum::color(SpectrumColor::White) * 0.1),
    );
    objects.push(big_light);

    // wavy bezier
    let start_x = 400.;
    let end_x = 800.;
    let y = 100.;

    let n_points = 30;
    let mut points = Vec::new();
    //points.push(DVec2::new(400., 100.));

    //let noise_loop = NoiseLoop::new(DVec2::ZERO, 1., 1., 10.);

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
    let bezier_obj = Object::new(Shape::Strip(bezier_strip), Material::Reflective);
    objects.push(bezier_obj);

    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(400., 0.), DVec2::new(800., 0.))),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );

    // ?
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(0., 250.), DVec2::new(100., 250.))),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(110., 250.), DVec2::new(200., 250.))),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
    );
    objects.push(wall);
    let wall = Object::new(
        Shape::Segment(Segment::new(DVec2::new(200., 250.), DVec2::new(200., 500.))),
        Material::diffuse(Spectrum::color(SpectrumColor::White)),
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
        Material::diffuse(Spectrum::color(SpectrumColor::Red)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 2. * d_phi) * r,
        center + DVec2::from_angle(angle + 1. * d_phi) * r,
        Material::diffuse(Spectrum::color(SpectrumColor::Green)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 3. * d_phi) * r,
        center + DVec2::from_angle(angle + 2. * d_phi) * r,
        Material::diffuse(Spectrum::color(SpectrumColor::Blue)),
    );
    objects.push(segment);

    let segment = Object::segment(
        center + DVec2::from_angle(angle + 4. * d_phi) * r,
        center + DVec2::from_angle(angle + 3. * d_phi) * r,
        Material::diffuse(Spectrum::color(SpectrumColor::Yellow)),
    );
    objects.push(segment);

    let dieletric = Object::from_points(
        &[
            DVec2::new(50., 300.),
            DVec2::new(200., 300.),
            DVec2::new(125., 350.),
        ],
        //Material::dieletric(1.5),
        Material::Dielectric {
            ior: Ior::Cauchy { a: 1.45, b: 0.0354 },
        },
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
