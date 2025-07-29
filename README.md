# Mitsuba0.6 OHD Renderer

### [Project Page](https://juhyeonkim95.github.io/project-pages/ohd_rendering/) | [Paper](https://dl.acm.org/doi/10.1145/3731150) | [Tutorial](ohd_tutorial/README.md)

![visualization](assets/teaser.png)

This repository is the official Mitsuba0.6 implementation of "A Monte Carlo Rendering Framework for Simulating Optical Heterodyne Detection" by 
[Juhyeon Kim](https://juhyeonkim.netlify.app/), 
[Craig Benko](http://craigbenko.com/), 
[Magnus Wrenninge](https://scholar.google.se/citations?user=aa85HZgAAAAJ&hl=en), 
[Ryusuke Villemin](https://www.linkedin.com/in/ryusuke-v-04649b1/), 
[Zeb Barber](https://www.linkedin.com/in/zeb-barber-a03020127/), 
[Wojciech Jarosz](https://cs.dartmouth.edu/~wjarosz/), 
[Adithya Pediredla](https://sites.google.com/view/adithyapediredla/)
(SIGGRAPH 2025, journal paper, 🏆 Honorable mention).

## Install
To compile, follow the original Mitsuba compliation guide at [here](https://github.com/mitsuba-renderer/mitsuba).

Instead of original config, please use double precision and mono-channel config (`config_double_single.py`).

## Parameter Explanation
New integrators are added in `src/integrators/fmcw` folder, for OHD rendering.
Two types of laser are typically used for OHD rendering, standard laser with constant frequency and swept-frequency laser with a chirp.
FMCW (frequency-modulated continuous-wave) use swept-frequency laser with up and down chirp.
Since FMCW is the most general case of OHD, we named integrator as `fmcw`.

### Integrator Types
* `fmcw_psd`: frequency-domain PSD (Power Spectral Density) rendering, corresponding to Algorithm1 in the paper. It outputs [H x W x 2M] where M is histogram bin number. 2M is for (up, down) chirp each.
* `fmcw_field`: time-domain rendering with field sampling and random phase addition, corresponding to Algorithm2 in the paper. It outputs [H x W x 4M] where M is FFT sample number. 4M is for (up, down) chirp x (real, complex) term each. 

### Detailed Parameters
* `wavelength` : Laser wavelength in nanometer (default : 1550)
* `f_0` : Laser frequency in GHz. If `wavelength` is also defined, we use `wavelength` first (default : 193,414)
* `T` : Sweep time for one chirp in microsecond (default : 10)
* `B` : Signal bandwidth in GHz (default : 1)
* `M` : Number of time stamps used for FFT or histogram bin (default : 4096)
* `use_collimated` : Ignore light setting and use collimated laser (default : false)
* `fov_error` : Additional parameter to change laser FOV. If set to $a$, instead of uniform sample from whole pixel, we sample centered square-area with $2a$ (default : 0.5).
* `use_random_phase` : Whether to use random phase addition for `fmcw_field`. (default : true)

## Usage
We included exhaustive examples in `ohd_tutorial` folder at [here](ohd_tutorial/README.md). 
You can simulate various experiments in the main paper including below interactive PSD visualization.
![interactive_fmcw_psd](ohd_tutorial/assets/interactive_fmcw_psd.gif)


## Citation
If you find this useful for your research, please consider to cite:
```
@article{kim2025ohd,
  author = {Kim, Juhyeon and Benko, Craig and Wrenninge, Magnus and Villemin, Ryusuke and Barber, Zeb and Jarosz, Wojciech and Pediredla, Adithya},
  title = {A Monte Carlo Rendering Framework for Simulating Optical Heterodyne Detection},
  year = {2025},
  issue_date = {August 2025},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {44},
  number = {4},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3731150},
  doi = {10.1145/3731150},
  journal = {ACM Trans. Graph.},
  month = jul,
  articleno = {56},
  numpages = {19}
}
```
