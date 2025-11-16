
rm -rf out/*.png
cargo r --release --bin 2d
ffmpeg.exe -y -framerate 30 -i ".\out\out_%04d.png" -c:v libx264 -vf fps=25 -pix_fmt yuv420p out.mp4
