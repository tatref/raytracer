use glam::{IVec2, Vec3};
use image::{ImageBuffer, Rgba};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

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

#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct PixelData {
    pub weight: f32,
    pub value: Vec3,
}

impl Default for PixelData {
    fn default() -> Self {
        PixelData {
            weight: 1.,
            value: Vec3::ZERO,
        }
    }
}

use std::ops::{Add, AddAssign, Div, Mul, Sub};
impl Add for PixelData {
    type Output = PixelData;
    fn add(self, rhs: Self) -> Self::Output {
        let mut pixel_data = PixelData::default();
        pixel_data.weight = self.weight + rhs.weight;
        pixel_data.value =
            (self.value * self.weight + rhs.value * rhs.weight) / (self.weight + rhs.weight);
        pixel_data
    }
}
impl AddAssign for PixelData {
    fn add_assign(&mut self, rhs: Self) {
        self.weight += rhs.weight;
        self.value =
            (self.value * self.weight + rhs.value * rhs.weight) / (self.weight + rhs.weight);
    }
}
impl Mul<f32> for PixelData {
    type Output = PixelData;
    fn mul(self, rhs: f32) -> Self::Output {
        let mut pixel_data = PixelData::default();
        pixel_data.value = self.value * rhs;
        pixel_data
    }
}
impl Div<f32> for PixelData {
    type Output = PixelData;
    fn div(self, rhs: f32) -> Self::Output {
        let mut pixel_data = PixelData::default();
        pixel_data.value = self.value / rhs;
        pixel_data
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct RawImage {
    pub width: i32,
    pub height: i32,
    pub data: Vec<PixelData>,
}
impl RawImage {
    pub fn new(width: i32, height: i32) -> Self {
        let data = vec![PixelData::default(); width as usize * height as usize];
        Self {
            width,
            height,
            data,
        }
    }

    fn save(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut f = std::fs::File::create(path)?;
        rmp_serde::encode::write(&mut f, self)?;

        Ok(())
    }

    pub fn pixel_to_idx(&self, pixel: IVec2) -> Result<usize, ()> {
        if pixel.x < 0 || pixel.x >= self.width || pixel.y < 0 || pixel.y >= self.height {
            return Err(());
        }

        Ok((self.width * pixel.y + pixel.x) as usize)
    }

    pub fn idx_to_pixel(&self, idx: usize) -> Result<IVec2, ()> {
        if idx > (self.width * self.height * 3) as usize {
            return Err(());
        }

        let pixel_index = idx;
        let y = pixel_index as i32 / self.width;
        let x = pixel_index as i32 % self.width;
        Ok(IVec2::new(x, y))
    }

    pub fn get(&self, pixel: IVec2) -> Result<PixelData, ()> {
        let idx = self.pixel_to_idx(pixel)?;
        Ok(self.data[idx])
    }

    pub fn draw_pixel(
        &mut self,
        pixel: IVec2,
        pixel_data: PixelData,
        blending: Blending,
    ) -> Result<(), ()> {
        let idx = self.pixel_to_idx(pixel)?;

        match blending {
            Blending::Add => {
                self.data[idx + 0] += pixel_data;
            }
            Blending::Replace => {
                self.data[idx + 0] = pixel_data;
            }
        }

        Ok(())
    }

    pub fn abs(&self) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (src, dst) in self.data.iter().zip(image.data.iter_mut()) {
            *dst = PixelData {
                weight: src.weight,
                value: src.value.abs(),
            }
        }

        image
    }

    pub fn map_pixel<T: Fn(PixelData) -> PixelData>(&mut self, f: T) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (src, dst) in self.data.iter().zip(image.data.iter_mut()) {
            *dst = f(*src);
        }

        image
    }

    pub fn convolution<const N: usize>(&self, kernel: [[f32; N]; N]) -> RawImage {
        let mut image = RawImage::new(self.width, self.height);

        for (x, y) in (0..(self.width - N as i32)).cartesian_product(0..(self.height - N as i32)) {
            let mut pixel_data = PixelData::default();

            for (i, j) in (0..(kernel.len() as i32)).cartesian_product(0..(kernel.len() as i32)) {
                let pixel = IVec2::new(x + i, y + j);
                pixel_data += self.get(pixel).unwrap_or(PixelData::default())
                    * kernel[i as usize][j as usize];
            }

            image
                .draw_pixel(IVec2::new(x, y), pixel_data, Blending::Replace)
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

        //let max = self.data.iter().cloned().reduce(f32::max).unwrap();
        let max = 0.;

        for (dst, src) in img.pixels_mut().zip(self.data.iter()) {
            match *tone_mapping_method {
                ToneMappingMethod::Gamma { gamma } => {
                    //dst = ((src / max).powf(gamma) * 255.) as u8;
                    unimplemented!();
                }
                ToneMappingMethod::Reinhard => {
                    dst[0] = (tone_mapping_reinhard(src.value.x) * 255.) as u8;
                    dst[1] = (tone_mapping_reinhard(src.value.y) * 255.) as u8;
                    dst[2] = (tone_mapping_reinhard(src.value.z) * 255.) as u8;
                }
            }

            dst[3] = 255;
        }

        img
    }
}

impl Add for RawImage {
    type Output = RawImage;
    fn add(self, rhs: Self) -> Self::Output {
        let mut sum = RawImage::new(self.width, self.height);

        for ((dest, src_self), src_rhs) in sum
            .data
            .iter_mut()
            .zip(self.data.iter())
            .zip(rhs.data.iter())
        {
            *dest = *src_self + *src_rhs;
        }

        sum
    }
}

impl Sub for RawImage {
    type Output = RawImage;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut sum = RawImage::new(self.width, self.height);

        for ((dest, src_self), src_rhs) in sum
            .data
            .iter_mut()
            .zip(self.data.iter())
            .zip(rhs.data.iter())
        {
            *dest = PixelData {
                value: src_self.value - src_rhs.value,
                weight: 1.,
            };
        }

        sum
    }
}
