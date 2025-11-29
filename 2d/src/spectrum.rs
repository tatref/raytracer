use glam::Vec3;
use palette::xyz::Xyz;

/// Data from http://www.cvrl.org/cmfs.htm, CIE 1931 2-deg, XYZ CMFs
/// 360 nm -> 830 nm += 5nm
const SPECTRUM_TO_XYZ_MAP: [[f32; 3]; 95] = [
    [0.000129900000, 0.000003917000, 0.000606100000],
    [0.000232100000, 0.000006965000, 0.001086000000],
    [0.000414900000, 0.000012390000, 0.001946000000],
    [0.000741600000, 0.000022020000, 0.003486000000],
    [0.001368000000, 0.000039000000, 0.006450001000],
    [0.002236000000, 0.000064000000, 0.010549990000],
    [0.004243000000, 0.000120000000, 0.020050010000],
    [0.007650000000, 0.000217000000, 0.036210000000],
    [0.014310000000, 0.000396000000, 0.067850010000],
    [0.023190000000, 0.000640000000, 0.110200000000],
    [0.043510000000, 0.001210000000, 0.207400000000],
    [0.077630000000, 0.002180000000, 0.371300000000],
    [0.134380000000, 0.004000000000, 0.645600000000],
    [0.214770000000, 0.007300000000, 1.039050100000],
    [0.283900000000, 0.011600000000, 1.385600000000],
    [0.328500000000, 0.016840000000, 1.622960000000],
    [0.348280000000, 0.023000000000, 1.747060000000],
    [0.348060000000, 0.029800000000, 1.782600000000],
    [0.336200000000, 0.038000000000, 1.772110000000],
    [0.318700000000, 0.048000000000, 1.744100000000],
    [0.290800000000, 0.060000000000, 1.669200000000],
    [0.251100000000, 0.073900000000, 1.528100000000],
    [0.195360000000, 0.090980000000, 1.287640000000],
    [0.142100000000, 0.112600000000, 1.041900000000],
    [0.095640000000, 0.139020000000, 0.812950100000],
    [0.057950010000, 0.169300000000, 0.616200000000],
    [0.032010000000, 0.208020000000, 0.465180000000],
    [0.014700000000, 0.258600000000, 0.353300000000],
    [0.004900000000, 0.323000000000, 0.272000000000],
    [0.002400000000, 0.407300000000, 0.212300000000],
    [0.009300000000, 0.503000000000, 0.158200000000],
    [0.029100000000, 0.608200000000, 0.111700000000],
    [0.063270000000, 0.710000000000, 0.078249990000],
    [0.109600000000, 0.793200000000, 0.057250010000],
    [0.165500000000, 0.862000000000, 0.042160000000],
    [0.225749900000, 0.914850100000, 0.029840000000],
    [0.290400000000, 0.954000000000, 0.020300000000],
    [0.359700000000, 0.980300000000, 0.013400000000],
    [0.433449900000, 0.994950100000, 0.008749999000],
    [0.512050100000, 1.000000000000, 0.005749999000],
    [0.594500000000, 0.995000000000, 0.003900000000],
    [0.678400000000, 0.978600000000, 0.002749999000],
    [0.762100000000, 0.952000000000, 0.002100000000],
    [0.842500000000, 0.915400000000, 0.001800000000],
    [0.916300000000, 0.870000000000, 0.001650001000],
    [0.978600000000, 0.816300000000, 0.001400000000],
    [1.026300000000, 0.757000000000, 0.001100000000],
    [1.056700000000, 0.694900000000, 0.001000000000],
    [1.062200000000, 0.631000000000, 0.000800000000],
    [1.045600000000, 0.566800000000, 0.000600000000],
    [1.002600000000, 0.503000000000, 0.000340000000],
    [0.938400000000, 0.441200000000, 0.000240000000],
    [0.854449900000, 0.381000000000, 0.000190000000],
    [0.751400000000, 0.321000000000, 0.000100000000],
    [0.642400000000, 0.265000000000, 0.000049999990],
    [0.541900000000, 0.217000000000, 0.000030000000],
    [0.447900000000, 0.175000000000, 0.000020000000],
    [0.360800000000, 0.138200000000, 0.000010000000],
    [0.283500000000, 0.107000000000, 0.000000000000],
    [0.218700000000, 0.081600000000, 0.000000000000],
    [0.164900000000, 0.061000000000, 0.000000000000],
    [0.121200000000, 0.044580000000, 0.000000000000],
    [0.087400000000, 0.032000000000, 0.000000000000],
    [0.063600000000, 0.023200000000, 0.000000000000],
    [0.046770000000, 0.017000000000, 0.000000000000],
    [0.032900000000, 0.011920000000, 0.000000000000],
    [0.022700000000, 0.008210000000, 0.000000000000],
    [0.015840000000, 0.005723000000, 0.000000000000],
    [0.011359160000, 0.004102000000, 0.000000000000],
    [0.008110916000, 0.002929000000, 0.000000000000],
    [0.005790346000, 0.002091000000, 0.000000000000],
    [0.004109457000, 0.001484000000, 0.000000000000],
    [0.002899327000, 0.001047000000, 0.000000000000],
    [0.002049190000, 0.000740000000, 0.000000000000],
    [0.001439971000, 0.000520000000, 0.000000000000],
    [0.000999949300, 0.000361100000, 0.000000000000],
    [0.000690078600, 0.000249200000, 0.000000000000],
    [0.000476021300, 0.000171900000, 0.000000000000],
    [0.000332301100, 0.000120000000, 0.000000000000],
    [0.000234826100, 0.000084800000, 0.000000000000],
    [0.000166150500, 0.000060000000, 0.000000000000],
    [0.000117413000, 0.000042400000, 0.000000000000],
    [0.000083075270, 0.000030000000, 0.000000000000],
    [0.000058706520, 0.000021200000, 0.000000000000],
    [0.000041509940, 0.000014990000, 0.000000000000],
    [0.000029353260, 0.000010600000, 0.000000000000],
    [0.000020673830, 0.000007465700, 0.000000000000],
    [0.000014559770, 0.000005257800, 0.000000000000],
    [0.000010253980, 0.000003702900, 0.000000000000],
    [0.000007221456, 0.000002607800, 0.000000000000],
    [0.000005085868, 0.000001836600, 0.000000000000],
    [0.000003581652, 0.000001293400, 0.000000000000],
    [0.000002522525, 0.000000910930, 0.000000000000],
    [0.000001776509, 0.000000641530, 0.000000000000],
    [0.000001251141, 0.000000451810, 0.000000000000],
];

pub const SPECTRUM_SAMPLES: usize = SPECTRUM_TO_XYZ_MAP.len();

/// ! Spectrum from 360 nm to 830 nm, increment by 5 nm
#[derive(Clone)]
pub struct Spectrum {
    pub data: [f32; SPECTRUM_SAMPLES],
}

impl Default for Spectrum {
    fn default() -> Self {
        let data = [0.; _];
        Self { data }
    }
}

#[derive(Clone, Copy, Debug)]
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
            .map(|(idx, power)| (360 + idx as u16 * 5, power))
    }

    pub fn color(color: SpectrumColor) -> Self {
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
        let mut x = 0.;
        let mut y = 0.;
        let mut z = 0.;

        for (intensity, map) in self.data.iter().zip(SPECTRUM_TO_XYZ_MAP.iter()) {
            x = x + *intensity * map[0];
            y = y + *intensity * map[1];
            z = z + *intensity * map[2];
        }

        Vec3::new(x, y, z)
    }
}

use std::ops::{Add, Div, DivAssign, Mul, MulAssign};
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
