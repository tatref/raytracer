#![allow(dead_code)]
#![allow(unused)]

use std::{
    io::Write,
    sync::mpsc::{self, Receiver, Sender, SyncSender},
    thread::{self, JoinHandle, Thread},
};

use eframe::egui::{self, CollapsingHeader, Image, TextureHandle, TextureOptions, Ui};
use enum2egui::GuiInspect;
use image::EncodableLayout;
use raytracer::{
    img::{RawImage, ToneMappingMethod},
    librt2d::*,
    spectrum::{Spectrum, SpectrumColor},
    worlds::*,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct CheckPoint {
    world: World,
    raw_image: RawImage,
}

impl CheckPoint {
    fn save_checkpoint(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::File::create(path)?;

        let serialized = rmp_serde::encode::to_vec(self)?;
        f.write_all(&serialized)?;

        Ok(())
    }

    fn load_checkpoint(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut f = std::fs::File::open(path)?;
        let val: CheckPoint = rmp_serde::decode::from_read(f)?;
        Ok(val)
    }
}

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size((1024., 768.)),
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
    CornellBoxAbsorption,
    CornellBoxSSS,
    Simple,
    Complex,
    Sample,
    Spectrum,
}

#[derive(Debug, PartialEq, Eq)]
struct ColorMetadata {
    color: SpectrumColor,
}

struct MyApp<'a> {
    load_world: WorldList,
    world: World,
    world_time: f64,
    texture_handle: Option<TextureHandle>,
    current_image: Option<(Image<'a>, RawImage)>,
    tx: Option<SyncSender<RenderCommand>>,
    rx: Option<Receiver<RenderProgress>>,
    render_thread: Option<JoinHandle<()>>,
    current_loop: usize,
    objects_metadata: Vec<ColorMetadata>,
    tone_mapping_param: f32,
}

impl<'a> Default for MyApp<'a> {
    fn default() -> Self {
        let width = 800;
        let height = 600;
        let spp = 100;
        let recursion_limit = 10;

        let denoiser = Denoiser {
            mask_size: 2,
            oversampling_factor: 1.,
            top: 0.001,
            passes: 0,
        };
        let render_params = RenderParams {
            height,
            spp,
            stop_condition: StopCondition::Endless,
            width,
            recursion_limit,
            lambda_samples: 2,
            denoiser: Some(denoiser),
            use_quadtree: false,
        };

        let world = spectrum_world(render_params, 0., 0);
        let objects_metadata = Self::refresh_objects_metadata(&world);

        Self {
            load_world: WorldList::Spectrum,
            world,
            world_time: 0.,
            texture_handle: None,
            current_image: None,
            tx: None,
            rx: None,
            render_thread: None,
            current_loop: 0,
            objects_metadata,
            tone_mapping_param: 1.,
        }
    }
}

impl<'a> MyApp<'a> {
    fn refresh_objects_metadata(world: &World) -> Vec<ColorMetadata> {
        let mut v = Vec::new();

        for obj in &world.objects {
            v.push(ColorMetadata {
                color: SpectrumColor::Black,
            });
        }

        v
    }

    fn render(&mut self) {
        println!("Starting render...");

        let (tx1, rx1) = mpsc::sync_channel(5);
        let (tx2, rx2) = mpsc::sync_channel(5);

        self.rx = Some(rx1);
        self.tx = Some(tx2);

        let mut world = self.world.clone();
        world.build_quadtree();
        let render_thread = thread::spawn(move || {
            world.endless_render(tx1, rx2);
        });
        self.render_thread = Some(render_thread);
    }
}

fn render_params_ui(ui: &mut Ui, render_params: &mut RenderParams) {
    ui.heading("Render params");
    ui.horizontal(|ui| {
        ui.label("spp");
        ui.add(egui::DragValue::new(&mut render_params.spp).speed(1));
    });

    render_params.stop_condition.ui_mut(ui);

    ui.checkbox(&mut render_params.use_quadtree, "quadtree");
    ui.horizontal(|ui| {
        ui.label("recursion limit");
        ui.add(egui::DragValue::new(&mut render_params.recursion_limit).speed(1));
    });
    ui.horizontal(|ui| {
        ui.label("lambda samples");
        ui.add(egui::DragValue::new(&mut render_params.lambda_samples).speed(1));
    });

    match &mut render_params.denoiser {
        Some(denoiser) => {
            ui.heading("Denoiser");
            ui.horizontal(|ui| {
                ui.label("passes");
                ui.add(egui::DragValue::new(&mut denoiser.passes).speed(1));
            });
            ui.horizontal(|ui| {
                ui.label("top");
                ui.add(egui::DragValue::new(&mut denoiser.top).speed(0.01));
            });
            ui.horizontal(|ui| {
                ui.label("oversampling factor");
                ui.add(egui::DragValue::new(&mut denoiser.oversampling_factor).speed(1));
            });
            ui.horizontal(|ui| {
                ui.label("mask size");
                ui.add(egui::DragValue::new(&mut denoiser.mask_size).speed(1));
            });
        }
        None => (),
    }
}

impl<'a> eframe::App for MyApp<'a> {
    //fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show_inside(ui, |ui| {
            egui::Panel::left("left_panel")
                .resizable(true)
                .default_size(150.0)
                .show_inside(ui, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.label("Time");
                        ui.add(
                            egui::DragValue::new(&mut self.world_time)
                                .speed(0.01)
                                .range(0..=1),
                        );

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
                                    WorldList::CornellBoxAbsorption,
                                    "Cornell Box Absorption",
                                );
                                ui.selectable_value(
                                    &mut self.load_world,
                                    WorldList::CornellBoxSSS,
                                    "Cornell Box SSS",
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
                                    WorldList::Spectrum,
                                    "Spectrum",
                                );
                            });

                        if ui.button("Load").clicked() {
                            let load_world_fn = match self.load_world {
                                WorldList::Spectrum => spectrum_world,
                                WorldList::CornellBox => cornell_box,
                                WorldList::CornellBoxAbsorption => cornell_box_absorption,
                                WorldList::CornellBoxSSS => cornell_box_sss,
                                WorldList::Simple => simple_world,
                                WorldList::Complex => complex_world,
                                WorldList::Sample => sample_world,
                            };

                            self.world =
                                load_world_fn(self.world.render_params.clone(), self.world_time, 0);

                            let objects_metadata = Self::refresh_objects_metadata(&self.world);
                        }

                        ui.separator();

                        render_params_ui(ui, &mut self.world.render_params);

                        if self.render_thread.is_none() {
                            if ui.button("Render").clicked() {
                                self.render();
                            }
                        } else {
                            if ui.button("Stop").clicked() {
                                self.tx.as_ref().unwrap().send(RenderCommand::StopRender);
                                self.render_thread.take().unwrap().join().unwrap();
                            }
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
                                Material::Emissive { emission } => {
                                    ui.label("Emission color");

                                    egui::ComboBox::from_id_salt(idx)
                                        .selected_text(format!(
                                            "{:?}",
                                            self.objects_metadata[idx].color
                                        ))
                                        .show_ui(ui, |ui| {
                                            // TODO
                                            ui.selectable_value(
                                                &mut self.objects_metadata[idx].color,
                                                SpectrumColor::Black,
                                                "Black",
                                            );
                                            ui.selectable_value(
                                                &mut self.objects_metadata[idx].color,
                                                SpectrumColor::Blue,
                                                "Blue",
                                            );
                                            ui.selectable_value(
                                                &mut self.objects_metadata[idx].color,
                                                SpectrumColor::White,
                                                "White",
                                            );
                                        });

                                    if ui.button("OK").clicked() {
                                        *emission = Spectrum::emission_from_color(
                                            self.objects_metadata[idx].color,
                                        );
                                    }
                                }
                                Material::Dielectric { ior, absorption } => {
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

            egui::CentralPanel::default().show_inside(ui, |ui| {
                ui.request_repaint_after_secs(0.2);

                ui.label("Tone mapping");
                ui.add(
                    egui::DragValue::new(&mut self.tone_mapping_param)
                        .speed(0.1)
                        .range(0..=1000),
                );

                self.rx.as_ref().map(|rx| {
                    if let Ok(render_progress) = rx.try_recv() {
                        self.current_loop = render_progress.loops;

                        let image = render_progress.raw_image.convert_to_image(
                            &ToneMappingMethod::Reinhard {
                                param: self.tone_mapping_param,
                            },
                        );

                        let image = egui::ColorImage::from_rgba_unmultiplied(
                            [image.width() as usize, image.height() as usize],
                            image.as_bytes(),
                        );
                        let image_size = image.size;

                        let handle =
                            ui.load_texture("image texture", image, TextureOptions::default());
                        self.texture_handle = Some(handle.clone());
                        let sized_image = egui::load::SizedTexture::new(
                            handle.id(),
                            egui::vec2(image_size[0] as f32, image_size[1] as f32),
                        );
                        let image = egui::Image::from_texture(sized_image);
                        self.current_image = Some((image, render_progress.raw_image));
                    };
                });

                match &self.current_image {
                    Some((image, raw_image)) => {
                        ui.label(format!("loop: {}", self.current_loop));
                        ui.add_sized(image.size().unwrap(), image.clone());

                        ui.separator();

                        if ui.button("Save out.png & checkpoint.dat").clicked() {
                            let image = raw_image.convert_to_image(&ToneMappingMethod::Reinhard {
                                param: self.tone_mapping_param,
                            });
                            let _ = image.save("out.png");

                            let checkpoint = CheckPoint {
                                raw_image: raw_image.clone(),
                                world: self.world.clone(),
                            };
                            checkpoint.save_checkpoint("checkpoint.dat").unwrap();
                        }
                    }
                    None => {}
                }
            });
        });
    }
}
