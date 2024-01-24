# Mitsuba0.6 FMCW Renderer Tutorial

There are three main files in `src` folder.

* `main_run_fmcw_rendering.py` : perform FMCW rendering for Cornell-Box scene. It generates image in `results` folder with size of $[H, W, 2M]$ where H, W are image size and $M$ is number of time stamps (used for FFT).
* `main_export_depth_velocity_map.py` : perform depth / velocity map reconstruction
* `main_plot_fft_interactive.py` : perform interactive FFT visualization

The result of `main_export_depth_velocity_map.py` looks like
![visualization](assets/depthmap.png)
![visualization](assets/velocitymap.png)



We will skip the details for the code, but interested reader could check out the detailed code.


```
cd fmcw_tutorial
python main_run_fmcw_rendering.py           # run FMCW rendering
python main_export_depth_velocity_map.py    # perform depth / velocity reconstruction
python main_plot_fft_interactive.py         # interactive FFT visualization
```