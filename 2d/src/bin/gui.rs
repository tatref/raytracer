#![allow(dead_code)]
#![allow(unused)]

use eframe::egui::{self, CollapsingHeader, Image, TextureHandle, TextureOptions, Ui};
use image::EncodableLayout;
use raytracer::{
    img::{RawImage, ToneMappingMethod},
    librt2d::*,
    worlds::*,
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

#[derive(Debug, PartialEq, Eq)]
enum WorldList {
    CornellBox,
    Simple,
    Complex,
    Sample,
    Colors,
}

struct MyApp<'a> {
    load_world: WorldList,
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

        let world = colors_world(render_params, 0., 0);

        Self {
            load_world: WorldList::Colors,
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
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.heading("Load World");
                        egui::ComboBox::from_label("Load world")
                            .selected_text(format!("{:?}", self.load_world))
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::CornellBox,
                                    "Cornell Box",
                                );
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::Complex,
                                    "Complex",
                                );
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::Sample,
                                    "Sample",
                                );
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::Simple,
                                    "Simple",
                                );
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::Colors,
                                    "Colors",
                                );
                            });

                        if ui.button("Load").clicked() {
                            self.world = match self.load_world {
                                WorldList::Colors => {
                                    colors_world(self.world.render_params.clone(), 0., 0)
                                }
                                WorldList::CornellBox => {
                                    cornell_box(self.world.render_params.clone(), 0., 0)
                                }
                                WorldList::Simple => {
                                    simple_world(self.world.render_params.clone(), 0., 0)
                                }
                                WorldList::Complex => {
                                    complex_world(self.world.render_params.clone(), 0., 0)
                                }
                                WorldList::Sample => {
                                    sample_world(self.world.render_params.clone(), 0., 0)
                                }
                            };
                        }

                        ui.separator();

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
                                                ui.add(
                                                    egui::DragValue::new(&mut circle.r).speed(1.0),
                                                );
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
                                _ => {
                                    ui.label("unknown object");
                                } //_ => unimplemented!(),
                            }

                            match &mut obj.mat {
                                Material::Emissive {
                                    emission: emission_color,
                                } => {
                                    ui.label("Emission color");
                                    //ui.horizontal(|ui| {
                                    //    ui.label("r");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut emission_color.x).speed(1),
                                    //    );
                                    //    ui.label("g");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut emission_color.y).speed(1),
                                    //    );
                                    //    ui.label("b");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut emission_color.z).speed(1),
                                    //    );
                                    //});
                                }
                                Material::Dielectric { ior } => {
                                    ui.label("ior");
                                    match ior {
                                        Ior::Simple(ior) => {
                                            ui.add(
                                                egui::DragValue::new(ior).speed(0.1).range(0..=5),
                                            );
                                        }
                                        _ => (),
                                    }
                                }
                                Material::Diffuse { absorption } => {
                                    ui.label("Absorption color");
                                    //ui.horizontal(|ui| {
                                    //    ui.label("r");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut absorption.x)
                                    //            .speed(0.1)
                                    //            .range(0..=1),
                                    //    );
                                    //    ui.label("g");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut absorption.y)
                                    //            .speed(0.1)
                                    //            .range(0..=1),
                                    //    );
                                    //    ui.label("b");
                                    //    ui.add(
                                    //        egui::DragValue::new(&mut absorption.z)
                                    //            .speed(0.1)
                                    //            .range(0..=1),
                                    //    );
                                    //});
                                }
                                _ => (),
                            }
                            ui.separator();
                        }
                    })
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
