# Mitsuba OHD Rendering Tutorial

## About
We provide several tutorials for OHD rendering.

### Speckle Simulation for Simple Case (in `tutorial_speckle_simple`)
This example provides a simple simulation on two peak scenario as shown in Fig.7 in the main paper.
`main.py` does
* Evaluate normalized covariance matrix
* Evaluate power spectrum (to show Alg1 and Alg2 converge to the same result)
* Evaluate intensity histogram at specific frequency
![speckle_simple](assets/image1.png)

### Doppler Simulation using Microsurface (in `tutorial_microsuface_doppler`, this uses GPU)
This example provides a microsurface simulation on Doppler effect (target vs spot velocity) which corresponds to Fig.10 in the main paper.
`main.py` does
* Evaluate PSD by generating microsurface (from Gaussian process) and actually moving it over each time stamp with field evaluation.
* Change condition of RMS height and correlation length (two parameters used for Gaussian surface generation) and compare result.
![microsuface_doppler](assets/image2.png)

### Rendering Simulation (in `tutorial_ohd_rendering`)
First, run `src/render_full.py` to render FMCW histogram. 
The output is [H, W, 2M] where 2M is for up and down chirp.
Then, run `src/visualizer_full.py` to visualize PSD for up and down chirp.
It runs interactively.
![](video.mp4)
