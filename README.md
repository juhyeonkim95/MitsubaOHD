# Mitsuba0.6 FMCW Renderer
This repository is the official Mitsuba0.6 implementation of "A Monte Carlo Rendering Framework for Frequency-Modulated Continuous-Wave LiDAR" (SIGGRAPH 2024 paper id 482)

## Install
To compile, follow the original Mitsuba compliation guide at [here](https://github.com/mitsuba-renderer/mitsuba).

Instead of config, please use double precision and mono-channel config (`config_double_single.py`).

## Parameter Explanation
New integrators are added in `src/integrators/fmcw` folder, for FMCW rendering.

### Integrator Type
* `fmcw`: Naive MC method (unbiased, $MN$ samples, output channel is $M$)
* `fmcw_static`: Correlated method for static scene (unbiased, $M$ samples, output channel is $M$)
* `fmcw_dynamic`: Correlated method + linear approximation for dynamic scene (biased, $2M$ samples, output channel is $2M$ with up,down chirp)

where $M$ is number of time-stamps used for FFT and $N$ is spp at each time stamp.


### Detailed Parameters
* `T` : Sweep time for one chirp in microsecond (default : 10)
* `B` : Signal bandwidth in GHz (default : 1)
* `M` : Number of time stamps used for FFT (default : 4096)
* `wavelength` : Laser wavelength in nanometer (default : 1550)
* `f_c` : Laser frequency in GHz. If `wavelength` is also defined, we use `wavelength` first (default : 77)
* `invert` : invert chirp-used for up/down chirp (default : false)
* `use_collimated` : Ignore light setting and use collimated laser (default : false)
* `fov_error` : Additional parameter to change laser FOV. If set to $a$, instead of uniform sample from whole pixel, we sample centered area with $2a$ (default : 0.5).
* `use_amplitude` : use sqrt for the output (default : false)
* `pdf_sqrt` : divide output with pdf sqrt (default : false)
* `spatial_correlation_mode`: method for path correlation (only used for `fmcw_dynamic`) (default : `ray_sampler`)


## Usage
```
mitsuba -L error config_example/example.xml
```
We also included exhaustive example configurations with result image.

## Citation
If you find this useful for your research, please consider to cite:
```
(TBA)
```
