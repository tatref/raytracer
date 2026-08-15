use raytracer::spectrum::{SPECTRUM_STEP_SIZE, SPECTRUM_TO_XYZ_MAP};

fn main() {
    let spectrum = SPECTRUM_TO_XYZ_MAP;

    let chunk_size = 10;

    for chunks in spectrum.chunks_exact(chunk_size) {
        let mut sum = [0., 0., 0.];
        for xyz in chunks {
            sum[0] += xyz[0];
            sum[1] += xyz[1];
            sum[2] += xyz[2];
        }

        let avg = [
            sum[0] / chunk_size as f32,
            sum[1] / chunk_size as f32,
            sum[2] / chunk_size as f32,
        ];

        println!("[{:.12}, {:.12}, {:.12}],", avg[0], avg[1], avg[2]);
    }

    let array_size = spectrum.len() / chunk_size;
    println!("array_size = {}", array_size);

    let step_size = SPECTRUM_STEP_SIZE * chunk_size as u16;
    println!("step_size = {}", step_size);
}
