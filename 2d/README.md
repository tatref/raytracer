
# Notes

Path tracing in production
https://jo.dreggn.org/path-tracing-in-production/2017/talk-jo/

Ray Tracing in One Weekend
https://raytracing.github.io/


# ffmpeg commands

```
rm -rf out/*.png
cargo r --release --bin 2d
ffmpeg -y -framerate 30 -i ".\out\out_%04d.png" -c:v libx264 -crf 15 -vf fps=30 -pix_fmt yuv420p out.mp4
```