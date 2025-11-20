
```
rm -rf out/*.png
cargo r --release --bin 2d
ffmpege -y -framerate 30 -i ".\out\out_%04d.png" -c:v libx264 -crf 15 -vf fps=30 -pix_fmt yuv420p out.mp4
```