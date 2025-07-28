from utils.gaussian_surface import *
from utils.visualizer import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from tqdm import tqdm, trange
import numpy as np
from scipy.fft import fft, ifft, fftfreq, fftshift
import matplotlib
matplotlib.use('Svg')
import os


def run_doppler_experiment(seed=1, correlation_length_coeff=1, rms_h_coeff=1, vel_coeff=1):
    """
        The function that simulates Doppler effect on Gaussian random surface. 
        The scene configuration is explained in the main paper Figure 10.
    """

    # object moving velocity
    vel = np.array([0, 0, vel_coeff])
    # light direction
    v1 = np.array([1,0,-1])
    # viewing direction
    v2 = np.array([-1,0,1])

    # unit vectors
    a1 = v1 / np.linalg.norm(v1)
    a2 = v2 / np.linalg.norm(v2)
    vm = (a2 - a1)

    # some constants
    c = 3e8
    wave_length = 1.55
    f_0 = c / wave_length

    # surface patch size (square shaped)
    patch_size = 400
    patch_N = 2048

    # sample area size (square shaped, should be much smaller than patch size)
    sample_size = 100
    sample_N = 1024 * int(max(rms_h_coeff, 1))
    
    # RMS & Correlation length
    rms_h = rms_h_coeff * wave_length
    correlation_length = correlation_length_coeff * wave_length
    
    # Generate Gaussian Random Surface
    surface = generate_gaussian_surface(rms_h, correlation_length, patch_size, patch_N, seed=seed)

    T = 10
    N_T = 3072
    sx = 0
    sy = 0
    Es = []

    # u2 : instantaenous velocity of the intersection point on the  macroscopic geometry along the ray
    # u1 : local coordinate change rate
    u2 = vel[2] / a1[2] * a1
    u1 = -vel + u2

    # move data to cuda
    surface_tensor = torch.tensor(surface, device='cuda').unsqueeze(0).unsqueeze(0)
    u2_tensor = torch.tensor(u2[None, None, :], device="cuda")
    vm_tensor = torch.tensor(vm[None, None, :], device="cuda")
    
    # Actually move microsurface over time
    dt = T / N_T
    sampled_surface = None
    for ti in range(N_T):
        rx = torch.rand((sample_N, sample_N), device="cuda") - 0.5
        ry = torch.rand((sample_N, sample_N), device="cuda") - 0.5

        sx += u1[0] * dt
        sy += u1[1] * dt

        # evalute heightmap values that belongs to sample region at current time stamp
        sampled_points = sample_heightmap_cuda(
            surface_tensor, patch_size, sample_size, 
            sx, sy, sample_N, rx, ry)
        
        # save sampled heightmap patch for visualization
        if ti == 0:
            sampled_surface = sampled_points[:,:,2].cpu().numpy()
            # plot_surface(sampled_points[:,:,2].cpu().numpy(), "1_%.2f.png"%rms_h_coeff, W=sample_size, vmin=-10, vmax=10)

        sampled_points += u2_tensor * dt * ti

        # evaluate field sum from sampled points
        phase = torch.sum(sampled_points * vm_tensor, dim=-1)
        phase = 2 * np.pi * phase / wave_length
        E = torch.sum(torch.exp(1j * phase)) * 1 / np.sqrt(sample_N * sample_N)
        Es.append(E.item())

    # perform FFT
    yf = np.asarray(Es)
    yf = fft(yf)
    yf = fftshift(yf)
    yf = np.abs(yf) ** 2 + 1e-10

    xf = fftfreq(N_T, T)
    xf = fftshift(xf) * N_T
    xf *= wave_length

    # eval GT target, spot velocity
    v_target = np.sum((a2 - a1) * vel, axis=0)
    v_spot = np.sum((a2 - a1) * u2, axis=0)

    return v_target, v_spot, xf, yf, surface, sampled_surface

def run_multiple_experiments(exp_number=1):
    # (1) change RMS Heights
    if exp_number == 1:
        cs = [10] # correlation length
        rs = [0.25, 0.5, 1.0, 2.0] # RMS height
        vs = [2] # velocity
    # (2) change Correlation Lengths
    elif exp_number == 2:
        cs = [40, 80] # correlation length
        rs = [2.0] # RMS height
        vs = [2] # velocity

    # number of surface realization used for ensemble-averaging
    N_realization = 64

    # create output folder
    output_folder = "result/exp_%d" % exp_number
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for c in cs:
        for r in rs:
            for v in vs:
                yf_avg = 0
                for i in trange(N_realization):
                    v_tg, v_spot, xf, yf, surface, sampled_surface = run_doppler_experiment(
                        seed=i + 1, correlation_length_coeff=c, rms_h_coeff=r, vel_coeff=v)
                    yf_avg += yf

                    # export height map
                    if i == 0:
                        plot_surface(sampled_surface, "%s/heightmap_c_%.2f_r_%.2f.png" % (output_folder, c, r), W=100, vmin=-10, vmax=10)

                # take ensemble average
                yf_avg /= N_realization
                yf_avg = 10 * np.log10(yf_avg) + 1e-2

                # plot PSD
                plt.plot(xf, yf_avg, linewidth=2)

                # plot gt target and spot velocities
                plt.axvline(x=v_tg, label="target", color='red', linestyle="--", linewidth=6, alpha=0.2)
                plt.axvline(x=v_spot, label="spot", color='blue', linestyle="--", linewidth=6, alpha=0.2)

                # plot format
                plt.xlim([-30, 30])
                plt.ylim([10, 90])
                plt.xlabel("velocity (prop to frequency)")
                plt.ylabel("power (dB)")
                
                plt.tight_layout()

                # export
                plt.savefig("%s/PSD_c_%.2f_r_%.2f.png" % (output_folder, c, r), dpi=600)
                # plt.savefig("%s/image_c_%.2f_r_%.2f.svg" % (output_folder, c, r), dpi=600)
                plt.close('all')

if __name__ == "__main__":
    run_multiple_experiments(1)
    run_multiple_experiments(2)
