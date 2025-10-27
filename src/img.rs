use glam::{UVec2, Vec3};
use image::{ImageBuffer, Rgba};

#[derive(Copy, Clone)]
pub enum ToneMappingMethod {
    Reinhard,
    Gamma { gamma: f32 },
}

#[derive(Copy, Clone)]
pub enum Blending {
    Add,
    Replace,
}

pub struct RawImage {
    pub width: u32,
    pub height: u32,
    data: Vec<f32>,
}
impl RawImage {
    pub fn new(width: u32, height: u32) -> Self {
        let data = vec![0f32; width as usize * height as usize * 3];
        Self {
            width,
            height,
            data,
        }
    }

    pub fn pixel_to_idx(&self, pixel: UVec2) -> usize {
        (self.width * pixel.y + pixel.x) as usize * 3
    }

    pub fn draw_pixel(&mut self, pixel: UVec2, color: Vec3, blending: &Blending) {
        let idx = self.pixel_to_idx(pixel);

        match blending {
            Blending::Add => {
                self.data[idx + 0] += color.x;
                self.data[idx + 1] += color.y;
                self.data[idx + 2] += color.z;
            }
            Blending::Replace => {
                self.data[idx + 0] = color.x;
                self.data[idx + 1] = color.y;
                self.data[idx + 2] = color.z;
            }
        }
    }

    pub fn convert_to_image(
        &self,
        tone_mapping_method: &ToneMappingMethod,
    ) -> ImageBuffer<Rgba<u8>, Vec<u8>> {
        let mut img: ImageBuffer<Rgba<u8>, Vec<u8>> = ImageBuffer::new(self.width, self.height);

        fn tone_mapping_reinhard(v: f32) -> f32 {
            v / (1. + v)
        }

        let max = self.data.iter().cloned().reduce(f32::max).unwrap();

        for (dst, src) in img.pixels_mut().zip(self.data.chunks_exact(3)) {
            match *tone_mapping_method {
                ToneMappingMethod::Gamma { gamma } => {
                    dst[0] = ((src[0] / max).powf(gamma) * 255.) as u8;
                    dst[1] = ((src[1] / max).powf(gamma) * 255.) as u8;
                    dst[2] = ((src[2] / max).powf(gamma) * 255.) as u8;
                }
                ToneMappingMethod::Reinhard => {
                    dst[0] = (tone_mapping_reinhard(src[0]) * 255.) as u8;
                    dst[1] = (tone_mapping_reinhard(src[1]) * 255.) as u8;
                    dst[2] = (tone_mapping_reinhard(src[2]) * 255.) as u8;
                }
            }

            dst[3] = 254;
        }

        img
    }
}
