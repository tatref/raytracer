use glam::{IVec2, Vec3};
use image::{ImageBuffer, Rgba};
use itertools::Itertools;

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

type Color = Vec3;

#[derive(Clone)]
pub struct RawImage {
    pub width: i32,
    pub height: i32,
    pub data: Vec<f32>,
}
impl RawImage {
    pub fn new(width: i32, height: i32) -> Self {
        let data = vec![0f32; width as usize * height as usize * 3];
        Self {
            width,
            height,
            data,
        }
    }

    pub fn pixel_to_idx(&self, pixel: IVec2) -> Result<usize, ()> {
        if pixel.x < 0 || pixel.x >= self.width || pixel.y < 0 || pixel.y >= self.height {
            return Err(());
        }

        Ok((self.width * pixel.y + pixel.x) as usize * 3)
    }

    pub fn idx_to_pixel(&self, idx: usize) -> Result<IVec2, ()> {
        if idx > (self.width * self.height * 3) as usize {
            return Err(());
        }

        let pixel_index = idx / 3;
        let y = pixel_index as i32 / self.width;
        let x = pixel_index as i32 % self.width;
        Ok(IVec2::new(x, y))
    }

    pub fn get(&self, pixel: IVec2) -> Result<Color, ()> {
        let idx = self.pixel_to_idx(pixel)?;
        Ok(Color::from_slice(&self.data[idx..(idx + 3)]))
    }

    pub fn draw_pixel(&mut self, pixel: IVec2, color: Vec3, blending: Blending) -> Result<(), ()> {
        let idx = self.pixel_to_idx(pixel)?;

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

        Ok(())
    }

    pub fn map<T: Fn(f32) -> f32>(&mut self, f: T) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (src, dst) in self.data.iter().zip(image.data.iter_mut()) {
            *dst = f(*src);
        }

        image
    }

    pub fn map_pixel<T: Fn(Vec3) -> Vec3>(&mut self, f: T) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (src, dst) in self
            .data
            .chunks_exact(3)
            .zip(image.data.chunks_exact_mut(3))
        {
            let p = Vec3::from_slice(src);
            let new_p = f(p);

            dst[0] = new_p.x;
            dst[1] = new_p.y;
            dst[2] = new_p.z;
        }

        image
    }

    pub fn convolution<const N: usize>(&self, kernel: [[f32; N]; N]) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (x, y) in (0..(self.width - N as i32)).cartesian_product(0..(self.height - N as i32)) {
            let mut color = Color::ZERO;

            for (i, j) in (0..(kernel.len() as i32)).cartesian_product(0..(kernel.len() as i32)) {
                let pixel = IVec2::new(x + i, y + j);
                color += self.get(pixel).unwrap_or(Color::ZERO) * kernel[i as usize][j as usize];
            }

            image
                .draw_pixel(IVec2::new(x, y), color, Blending::Replace)
                .unwrap();
        }

        image
    }

    pub fn convert_to_image(
        &self,
        tone_mapping_method: &ToneMappingMethod,
    ) -> ImageBuffer<Rgba<u8>, Vec<u8>> {
        let mut img: ImageBuffer<Rgba<u8>, Vec<u8>> =
            ImageBuffer::new(self.width as u32, self.height as u32);

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
