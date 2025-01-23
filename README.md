# Mitsuba0.6 FMCW Renderer
This repository is the official Mitsuba0.6 implementation of "A Monte Carlo Rendering Framework for Simulating Optical Heterodyne Detection" (SIGGRAPH 2025 paper id 247)
This repository serves FMCW simulation.

## Install
To compile, follow the original Mitsuba compliation guide at [here](https://github.com/mitsuba-renderer/mitsuba).

Instead of config, please use double precision and mono-channel config (`config_double_single.py`).

## Parameter Explanation
New integrators are added in `src/integrators/fmcw` folder, for FMCW rendering.

### Integrator Type
* `fmcw_instant`: time-domain rendering. It outputs [H x W x M] where M is FFT sample number.
* `fmcw_instant_frequency`: frequency-domain rendering. It outputs [H x W x M] where M is histogram bin number.


### Detailed Parameters
* `T` : Sweep time for one chirp in microsecond (default : 10)
* `B` : Signal bandwidth in GHz (default : 1)
* `M` : Number of time stamps used for FFT (default : 4096)
* `wavelength` : Laser wavelength in nanometer (default : 1550)
* `f_c` : Laser frequency in GHz. If `wavelength` is also defined, we use `wavelength` first (default : 77)
* `use_collimated` : Ignore light setting and use collimated laser (default : false)
* `fov_error` : Additional parameter to change laser FOV. If set to $a$, instead of uniform sample from whole pixel, we sample centered area with $2a$ (default : 0.5).


## Usage
```
(TBA)
```


## Citation
If you find this useful for your research, please consider to cite:
```
(TBA)
```
