#![allow(dead_code)]
#![allow(unused)]

use std::cell;

use eframe::egui::{self, CollapsingHeader, Image, TextureHandle, TextureOptions, Ui};
use egui_extras::{Column, TableBuilder};
use image::Rgba;
use palette::{IntoColor, Srgb};
use raytracer::spectrum::{SPECTRUM_SAMPLES, Spectrum, SpectrumColor};

fn to_srgb(spectrum: &Spectrum) -> Srgb<u8> {
    let xyz = spectrum.to_xyz();
    let srgb: Srgb = xyz.into_color();
    let srgb: Srgb<u8> = srgb.into_format();

    srgb
}

fn full_spectrum() {
    let n_colors = SpectrumColor::iter_colors().len() as u32;
    let cell_size = 10;

    let mut img: image::ImageBuffer<Rgba<u8>, Vec<u8>> =
        image::ImageBuffer::new(n_colors * cell_size, n_colors * cell_size);

    for (idx_a, color_a) in SpectrumColor::iter_colors().iter().enumerate() {
        let spectrum_a = Spectrum::color(*color_a);
        let test = to_srgb(&spectrum_a);
        println!("{:?}: {:?}", color_a, test);

        for (idx_b, color_b) in SpectrumColor::iter_colors().iter().enumerate() {
            let spectrum_b = Spectrum::color(*color_b);

            let spectrum = (spectrum_a.clone() + spectrum_b) * 100.;
            let srgb = to_srgb(&spectrum);

            for i in 0..cell_size {
                for j in 0..cell_size {
                    let px = img.get_pixel_mut(
                        idx_a as u32 * cell_size + i as u32,
                        idx_b as u32 * cell_size + j as u32,
                    );

                    px.0[0] = srgb.red;
                    px.0[1] = srgb.green;
                    px.0[2] = srgb.blue;
                    px.0[3] = 255;
                }
            }
        }
    }

    img.save("out.png").unwrap();
}

fn main() {
    full_spectrum();
}
