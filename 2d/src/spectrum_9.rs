use glam::Vec3;
use palette::{IntoColor, Srgb, xyz::Xyz};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Serialize, Serializer};

/// Data from http://www.cvrl.org/cmfs.htm, CIE 1931 2-deg, XYZ CMFs
/// 360 nm -> 830 nm += 5nm
pub const SPECTRUM_TO_XYZ_MAP: [[f32; 3]; 9] = [
    [0.005451550242, 0.000152129185, 0.025843411684],
    [0.243393018842, 0.018193000928, 1.231777906418],
    [0.108696006238, 0.184272006154, 0.825877070427],
    [0.219811990857, 0.832050025463, 0.052560001612],
    [0.896320044994, 0.817700028419, 0.001710000099],
    [0.604205012321, 0.263000011444, 0.000097999990],
    [0.057478010654, 0.021066403016, 0.000000000000],
    [0.001902146847, 0.000686900050, 0.000000000000],
    [0.000054891756, 0.000019822419, 0.000000000000],
];

pub const SPECTRUM_SAMPLES: usize = SPECTRUM_TO_XYZ_MAP.len();
pub const SPECTRUM_LAMBDA_MIN: u16 = 360;
pub const SPECTRUM_LAMBDA_MAX: u16 = 830;
pub const SPECTRUM_STEP_SIZE: u16 = 50;

/// ! Spectrum from 360 nm to 830 nm, increment by 5 nm
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Spectrum {
    pub data: [f32; SPECTRUM_SAMPLES],
}

impl<'de> Deserialize<'de> for Spectrum {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        todo!()
    }
}
impl Serialize for Spectrum {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(self.data.len()))?;
        for e in &self.data {
            seq.serialize_element(e)?;
        }
        seq.end()
    }
}

impl Default for Spectrum {
    fn default() -> Self {
        let data = [0.; _];
        Self { data }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectrumColor {
    Black,
    Violet,
    Blue,
    Cyan,
    Green,
    Yellow,
    Orange,
    Red,
    White,
}

impl SpectrumColor {
    pub fn iter_colors() -> Vec<SpectrumColor> {
        vec![
            SpectrumColor::Black,
            SpectrumColor::Violet,
            SpectrumColor::Blue,
            SpectrumColor::Cyan,
            SpectrumColor::Green,
            SpectrumColor::Yellow,
            SpectrumColor::Orange,
            SpectrumColor::Red,
            SpectrumColor::White,
        ]
    }
}

impl Spectrum {
    /// Construct a new `Spectrum` from an array
    /// Fails if some value of data is less than zero
    pub fn new(data: [f32; SPECTRUM_SAMPLES]) -> Result<Self, ()> {
        if data.iter().any(|&intensity| intensity < 0.) {
            Err(())
        } else {
            Ok(Spectrum { data: data })
        }
    }

    pub fn rand_lambda() -> (usize, f64) {
        let i = rand::random_range(0..SPECTRUM_SAMPLES);
        return (
            i,
            i as f64 * SPECTRUM_STEP_SIZE as f64 + SPECTRUM_LAMBDA_MIN as f64,
        );
    }

    pub fn from(data: &[f32]) -> Result<Self, ()> {
        if data.len() != SPECTRUM_SAMPLES {
            return Err(());
        }

        let mut spectrum = Spectrum::default();
        for (from, to) in data.iter().zip(spectrum.data.iter_mut()) {
            *to = *from;
        }

        Ok(spectrum)
    }

    /// (lambda (nm), power)
    pub fn iter_lambda(&mut self) -> impl Iterator<Item = (u16, &mut f32)> {
        self.data
            .iter_mut()
            .enumerate()
            .map(|(idx, power)| (SPECTRUM_LAMBDA_MIN + idx as u16 * SPECTRUM_STEP_SIZE, power))
    }

    pub fn absorption_from_color(color: SpectrumColor) -> Self {
        Self::emission_from_color(SpectrumColor::White) - Self::emission_from_color(color)
    }

    pub fn emission_from_color(color: SpectrumColor) -> Self {
        let lambda = match color {
            SpectrumColor::Black => return Spectrum::default(),
            SpectrumColor::Violet => (380, 450),
            SpectrumColor::Blue => (450, 485),
            SpectrumColor::Cyan => (485, 500),
            SpectrumColor::Green => (500, 565),
            SpectrumColor::Yellow => (565, 590),
            SpectrumColor::Orange => (590, 625),
            SpectrumColor::Red => (625, 750),
            SpectrumColor::White => (360, 830),
        };

        let mut spectrum = Spectrum::default();
        spectrum.iter_lambda().for_each(|(l, power)| {
            if l >= lambda.0 && l <= lambda.1 {
                //*power = 10. / (lambda.1 - lambda.0) as f32;
                *power = 1.;
            }
        });
        spectrum
    }

    /// Converts a `Spectrum` to `Xyz` tristimulis values
    pub fn to_xyz(&self) -> Xyz<palette::white_point::D65, f32> {
        let mut x = 0.;
        let mut y = 0.;
        let mut z = 0.;

        for (intensity, map) in self.data.iter().zip(SPECTRUM_TO_XYZ_MAP.iter()) {
            x = x + *intensity * map[0];
            y = y + *intensity * map[1];
            z = z + *intensity * map[2];
        }

        Xyz::new(x, y, z)
    }

    pub fn to_dvec3(&self) -> Vec3 {
        fn to_srgb(spectrum: &Spectrum) -> Srgb<f32> {
            let xyz = spectrum.to_xyz();
            let srgb: Srgb = xyz.into_color();
            //let srgb: Srgb<u8> = srgb.into_format();

            srgb
        }

        //let mut x = 0.;
        //let mut y = 0.;
        //let mut z = 0.;

        //for (intensity, map) in self.data.iter().zip(SPECTRUM_TO_XYZ_MAP.iter()) {
        //    x = x + *intensity * map[0];
        //    y = y + *intensity * map[1];
        //    z = z + *intensity * map[2];
        //}

        let srgb = to_srgb(self);

        Vec3::new(srgb.red, srgb.green, srgb.blue)
    }
}

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub};
impl Add for Spectrum {
    type Output = Spectrum;
    fn add(self, rhs: Self) -> Self::Output {
        let mut spectrum = Spectrum::default();
        for ((dest, a), b) in spectrum
            .data
            .iter_mut()
            .zip(self.data.iter())
            .zip(rhs.data.iter())
        {
            *dest = a + b;
        }
        spectrum
    }
}
impl Sub for Spectrum {
    type Output = Spectrum;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut spectrum = Spectrum::default();
        for ((dest, a), b) in spectrum
            .data
            .iter_mut()
            .zip(self.data.iter())
            .zip(rhs.data.iter())
        {
            *dest = a - b;
        }
        spectrum
    }
}
impl AddAssign for Spectrum {
    fn add_assign(&mut self, rhs: Self) {
        for (dest, a) in self.data.iter_mut().zip(rhs.data.iter()) {
            *dest += a;
        }
    }
}

impl Mul<f32> for Spectrum {
    type Output = Spectrum;
    fn mul(self, rhs: f32) -> Self::Output {
        let mut spectrum = Spectrum::default();
        for (dest, a) in spectrum.data.iter_mut().zip(self.data.iter()) {
            *dest = a * rhs;
        }
        spectrum
    }
}

impl Mul<Spectrum> for Spectrum {
    type Output = Spectrum;
    fn mul(self, rhs: Spectrum) -> Self::Output {
        let mut spectrum = Spectrum::default();
        for ((dest, a), b) in spectrum
            .data
            .iter_mut()
            .zip(self.data.iter())
            .zip(rhs.data.iter())
        {
            *dest = a * b;
        }
        spectrum
    }
}

impl MulAssign<f32> for Spectrum {
    fn mul_assign(&mut self, rhs: f32) {
        for dest in self.data.iter_mut() {
            *dest = *dest * rhs;
        }
    }
}

impl Div<f32> for Spectrum {
    type Output = Spectrum;
    fn div(self, rhs: f32) -> Self::Output {
        let mut spectrum = Spectrum::default();
        for (dest, a) in spectrum.data.iter_mut().zip(self.data.iter()) {
            *dest = a / rhs;
        }
        spectrum
    }
}

impl DivAssign<f32> for Spectrum {
    fn div_assign(&mut self, rhs: f32) {
        for dest in self.data.iter_mut() {
            *dest = *dest / rhs;
        }
    }
}
