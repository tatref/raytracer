use std::clone;

use eframe::egui::{self, CollapsingHeader, Image, TextureHandle, TextureOptions, Ui};
use glam::DVec2;
use image::{EncodableLayout, ImageBuffer, Rgba};
use raytracer::{
    Color,
    img::{RawImage, ToneMappingMethod},
    librt2d::*,
};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_maximized(true),
        ..Default::default()
    };
    eframe::run_native(
        "My egui App",
        options,
        Box::new(|cc| {
            // This gives us image support:

            egui_extras::install_image_loaders(&cc.egui_ctx);

            Ok(Box::<MyApp>::default())
        }),
    )
}

#[derive(Debug, PartialEq, Eq)]
enum ShapeKind {
    Circle,
    Segment,
    Strip,
}

#[derive(Debug, PartialEq, Eq)]
enum MaterialKind {
    Diffuse,
    Dieletric,
    Reflective,
}

struct MyApp<'a> {
    world: World,
    texture_handle: Option<TextureHandle>,
    current_image: Option<(Image<'a>, RawImage)>,
}

impl<'a> Default for MyApp<'a> {
    fn default() -> Self {
        let width = 800;
        let height = 600;
        let spp = 50;
        let recursion_limit = 10;

        let render_params = RenderParams {
            height,
            spp,
            width,
            recursion_limit,
            //denoiser,
            denoiser: None,
        };

        fn cornell_box(render_params: RenderParams, _t: f64, _idx: u64) -> World {
            let mut objects = Vec::new();

            let center = DVec2::new(
                render_params.width as f64 / 2.,
                render_params.height as f64 / 2.,
            );
            let light = Object::new(
                Shape::Circle(Circle::new(center, 5.)),
                Material::emissive_at(50., Color::ONE * 20.),
            );
            objects.push(light);

            let back = Object::segment(
                DVec2::new(200., 100.),
                DVec2::new(600., 100.),
                Material::diffuse(Color::ONE),
            );
            objects.push(back);

            let right = Object::segment(
                DVec2::new(600., 100.),
                DVec2::new(600., 500.),
                Material::diffuse(Color::Y),
            );
            objects.push(right);

            let left = Object::segment(
                DVec2::new(200., 500.),
                DVec2::new(200., 100.),
                Material::diffuse(Color::X),
            );
            objects.push(left);

            let sphere = Object::new(
                Shape::Circle(Circle::new(DVec2::new(300., 200.), 50.)),
                Material::diffuse(Color::ONE),
            );
            objects.push(sphere);

            let sphere = Object::new(
                Shape::Circle(Circle::new(DVec2::new(500., 400.), 50.)),
                Material::dieletric(0.7),
            );
            objects.push(sphere);

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

            let center = DVec2::new(500., 250.);
            let light = Object::new(
                Shape::Circle(Circle::new(center, 1.)),
                Material::emissive_at(50., Color::ONE * 20.),
            );
            objects.push(light);

            let segment = Object::segment(
                DVec2::new(300., 200.),
                DVec2::new(500., 200.),
                Material::diffuse(Color::ONE),
            );
            objects.push(segment);

            let world = World::new(objects, render_params);

            world
        }

        //let world = simple_world(render_params, 0., 0);
        let world = cornell_box(render_params, 0., 0);

        Self {
            world,
            texture_handle: None,
            current_image: None,
        }
    }
}

impl<'a> MyApp<'a> {
    fn render(&mut self, ctx: &egui::Context) {
        let render_params = RenderParams {
            height: 600,
            spp: 20,
            width: 800,
            recursion_limit: 10,
            //denoiser,
            denoiser: None,
        };
        let t = 0.;
        let chrono = std::time::Instant::now();
        let mut raw_image = self.world.render();
        let elapsed = chrono.elapsed();
        dbg!(elapsed);
        //annotate(&mut raw_image, &self.world, self.p);
        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);

        let image = egui::ColorImage::from_rgba_unmultiplied(
            [image.width() as usize, image.height() as usize],
            image.as_bytes(),
        );
        let image_size = image.size;

        let handle = ctx.load_texture("image texture", image, TextureOptions::default());
        self.texture_handle = Some(handle.clone());
        let sized_image = egui::load::SizedTexture::new(
            handle.id(),
            egui::vec2(image_size[0] as f32, image_size[1] as f32),
        );
        let image = egui::Image::from_texture(sized_image);
        self.current_image = Some((image, raw_image));
    }
}

fn render_params_ui(ui: &mut Ui, render_params: &mut RenderParams) {
    ui.heading("Render params");
    ui.horizontal(|ui| {
        ui.label("spp");
        ui.add(egui::DragValue::new(&mut render_params.spp).speed(1));
    });
    ui.horizontal(|ui| {
        ui.label("recursion limit");
        ui.add(egui::DragValue::new(&mut render_params.recursion_limit).speed(1));
    });
}

impl<'a> eframe::App for MyApp<'a> {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::SidePanel::left("left_panel")
                .resizable(true)
                .default_width(150.0)
                .show(ctx, |ui| {
                    render_params_ui(ui, &mut self.world.render_params);

                    if ui.button("Render").clicked() {
                        self.render(ctx);
                    }

                    ui.separator();

                    ui.heading("Objects");

                    for (idx, obj) in &mut self.world.objects.iter_mut().enumerate() {
                        match &mut obj.shape {
                            Shape::Circle(circle) => {
                                CollapsingHeader::new("Circle")
                                    .id_salt(idx)
                                    .default_open(false)
                                    .show(ui, |ui| {
                                        ui.horizontal(|ui| {
                                            ui.label("r");
                                            ui.add(egui::DragValue::new(&mut circle.r).speed(1.0));
                                        });
                                        ui.horizontal(|ui| {
                                            ui.label("x");
                                            ui.add(
                                                egui::DragValue::new(&mut circle.center.x)
                                                    .speed(1.0),
                                            );
                                            ui.label("y");
                                            ui.add(
                                                egui::DragValue::new(&mut circle.center.y)
                                                    .speed(1.0),
                                            );
                                        });
                                    });
                            }
                            Shape::Segment(segment) => {
                                ui.label("segment");
                            }
                            _ => unimplemented!(),
                        }

                        match &mut obj.mat {
                            Material::Emissive { emission_color, d } => {
                                ui.label("Emission color");
                                ui.horizontal(|ui| {
                                    ui.label("r");
                                    ui.add(egui::DragValue::new(&mut emission_color.x).speed(1));
                                    ui.label("g");
                                    ui.add(egui::DragValue::new(&mut emission_color.y).speed(1));
                                    ui.label("b");
                                    ui.add(egui::DragValue::new(&mut emission_color.z).speed(1));
                                });
                            }
                            Material::Dielectric { ior } => {
                                ui.label("ior");
                                ui.add(egui::DragValue::new(ior).speed(0.1).range(0..=5));
                            }
                            Material::Diffuse { absorption } => {
                                ui.label("Absorption color");
                                ui.horizontal(|ui| {
                                    ui.label("r");
                                    ui.add(
                                        egui::DragValue::new(&mut absorption.x)
                                            .speed(0.1)
                                            .range(0..=1),
                                    );
                                    ui.label("g");
                                    ui.add(
                                        egui::DragValue::new(&mut absorption.y)
                                            .speed(0.1)
                                            .range(0..=1),
                                    );
                                    ui.label("b");
                                    ui.add(
                                        egui::DragValue::new(&mut absorption.z)
                                            .speed(0.1)
                                            .range(0..=1),
                                    );
                                });
                            }
                            _ => (),
                        }
                        ui.separator();
                    }
                });

            //egui::SidePanel::right("right panel")
            //    .resizable(true)
            //    .default_width(900.)
            //    .show_inside();

            egui::CentralPanel::default().show(ctx, |ui| match &self.current_image {
                Some((image, raw_image)) => {
                    ui.add_sized(image.size().unwrap(), image.clone());

                    ui.separator();

                    if ui.button("Save out.png").clicked() {
                        let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard);
                        let _ = image.save("out.png");
                    }
                }
                None => {}
            });
        });
    }
}
