# Mitsuba OHD Rendering Tutorial

## About
We provide several tutorials for OHD rendering.
(still working on adding more tutorials)

## Non-rendering tutorial (this does not need Mitsuba)
### Speckle Simulation for Simple Case (in `tutorial_speckle_simple`)
This example provides a simple simulation on two peak scenario as shown in Fig.7 in the main paper.
`main.py` does
* Evaluate normalized covariance matrix
* Evaluate power spectrum (to show Alg1 and Alg2 converge to the same result)
* Evaluate intensity histogram at specific frequency
![speckle_simple](assets/speckle_simple.png)

### Doppler Simulation using Microsurface (in `tutorial_microsuface_doppler`, this uses GPU)
This example provides a microsurface simulation on Doppler effect (target vs spot velocity) which corresponds to Fig.10 in the main paper.
`main.py` does
* Evaluate PSD by generating microsurface (from Gaussian process) and actually moving it over each time stamp with field evaluation.
* Change condition of RMS height and correlation length (two parameters used for Gaussian surface generation) and compare result.
![microsuface_doppler](assets/microsuface_doppler.png)


## Rendering tutorial (this needs Mitsuba, in `tutorial_ohd_rendering`)
### Interactive histogram (using Alg1)
Run `src/render_full.py` to render FMCW PSD histogram (Alg1). 
The output is [H, W, 2M] where 2M is for up and down chirp.
Then, run `src/visualizer_full.py` to visualize PSD for up and down chirp.
It runs interactively.
![interactive_fmcw_psd](assets/interactive_fmcw_psd.gif)

### Single point comparison 
Run `src/render_single.py` to get FMCW PSD histogram (Alg1) and multiple field evaluations (Alg2). 
For PSD evaluation, output is [1, 1, 2M] where 2M is for up and down chirp.
For field evaluation, output is [1, 1, 4M] where 2M is for up and down chirp with complex value.
Then, run `src/visualizer_single.py` to visualize convergence of each algorithm for up and down chirp which corresponds to Fig.11 in the main paper. 
![rendering_algorithm_convergence](assets/rendering_algorithm_convergence.png)