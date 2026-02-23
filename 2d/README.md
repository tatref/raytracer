# What is it?
A 2D spectral raytracer, nearly functional



# Notes

Path tracing in production
https://jo.dreggn.org/path-tracing-in-production/2017/talk-jo/

Ray Tracing in One Weekend
https://raytracing.github.io/

Spectral Raytracing
https://graphics.cg.uni-saarland.de/courses/ris-2021/slides/Spectral%20Raytracing.pdf


# Gallery

Gui
![gui](gui.png)


Sample scene that can be animated. Featuring diffraction, bezier curves, mirrors...
![sample_scene](sample_scene.png)

![animated](animated.mp4)


Denoiser: computing difference between global image and current image
![densoise_laplace.png](denoise_laplace.png)


Denoiser: mask to select pixels to be recomputed
![denoise_mask.png](denoise_mask.png)


Denoiser: heatmap, the pixels with the most samples per pixel
![denoise_heatmap.png](denoise_heatmap.png)


Visible spectrum. Top: colored emission light. Bottom: white light + colored absorption box
![Visible spectrum](visible_spectrum2.png)


 Dispersion test. Left [Cachy's equation](https://en.wikipedia.org/wiki/Cauchy%27s_equation) with dispersion. Right: Diffraction without dispersion.
 ![diffraction_dispersion.png](diffration_dispersion.png)


Test scene with 200 spheres emission, absorption, diffraction. Probably faster some day with space partitioning
![stress_test](stress_test.png)


