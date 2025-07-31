import numpy as np
import matplotlib.pyplot as plt
import glob
from utils.exp_configs_utils import *
import matplotlib
import os
matplotlib.use("svg")

scene_name = "cornell-box-floor-specular"

# --- Configuration ---
base_name = "../results/%s/single/fmcwfield_single_depth_4/" % scene_name
file_pattern = base_name + "iter_*.npy"
file_list = sorted(glob.glob(file_pattern))
N = len(file_list)  # number of iteration files

if N == 0:
    raise ValueError("No matching files found!")

# Time and radar parameters
T = 5e-6
B = 1e9
c = 3e8

spp = 256
scene_scale = 10

# Load one file to determine M
sample_data = np.load(file_list[0]).squeeze()
M = sample_data.shape[0] // 4

# Distance axis for simulated FFT
freq = np.fft.fftfreq(M, T)
distance = np.fft.fftshift(freq) * M
distance = distance * c * T / (2 * B) * 2  # scale to distance

# FFT power for up/down chirp
up_powers = []
down_powers = []

for fname in file_list:
    raw = np.load(fname).squeeze()
    raw *= np.sqrt(spp)

    up = raw[0:M] + 1j * raw[2*M:3*M]
    down = raw[M:2*M] + 1j * raw[3*M:4*M]

    up_fft = np.fft.fftshift(np.fft.fft(up))
    down_fft = np.fft.fftshift(np.fft.fft(down))

    up_powers.append(np.abs(up_fft)**2 / (M))
    down_powers.append(np.abs(down_fft)**2 / (M))

up_powers = np.array(up_powers)
down_powers = np.array(down_powers)

# Average over all iterations
avg_up = up_powers.mean(axis=0)
avg_down = down_powers.mean(axis=0)
C_up = np.sum(avg_up) / M
C_down = np.sum(avg_up) / M
avg_up /= C_up
avg_down /= C_down


def to_log_scale(x, eps=1e-12):
    return 10 * np.log10(x + eps)

# Log-scale (dB)
log_up_sample = to_log_scale(up_powers[0] / C_up)
log_up_avg = to_log_scale(avg_up)
log_down_sample = to_log_scale(down_powers[0] / C_down)
log_down_avg = to_log_scale(avg_down)

# Flip down-chirp data for frequency reversal
log_down_sample_flipped = log_down_sample[::-1]
log_down_avg_flipped = log_down_avg[::-1]

# --- Load additional PSD data ---
psd_data = np.load("../results/%s/single/fmcwpsd_single_depth_4.npy" % scene_name).squeeze()  # shape: (2 × bin_size,)
bin_size = psd_data.shape[0] // 2
psd_data /= np.max(psd_data)

# Generate x-axis for PSD: (0, 400)
xmin, xmax = get_transient_min_max_range(scene_name)
psd_distance = np.linspace(xmin * scene_scale, xmax * scene_scale, bin_size)


# Extract up/down chirp PSD and flip down for frequency alignment
psd_up = psd_data[:bin_size]
psd_down = psd_data[bin_size:]
C_psd_up = np.sum(psd_up) / bin_size
C_psd_down = np.sum(psd_down) / bin_size

# normalize
psd_up /= C_psd_up
psd_down /= C_psd_down

# make distribution densor to reduce error with sinc
psd_distance_dense = np.linspace(xmin * scene_scale, xmax * scene_scale, bin_size * 4)
psd_up_dense = np.interp(psd_distance_dense, psd_distance, psd_up)
psd_down_dense = np.interp(psd_distance_dense, psd_distance, psd_down)

# perform convolution
sinc_kernel = np.sinc((distance[:, None] - psd_distance_dense[None, :]) / (distance[0] - distance[1])) ** 2
sinc_kernel /= np.sum(sinc_kernel, axis=1, keepdims=True)
psd_up_convolved = sinc_kernel @ psd_up_dense
psd_down_convolved = sinc_kernel @ psd_down_dense

# to log scale
psd_up = to_log_scale(psd_up)
psd_down = to_log_scale(psd_down)
psd_up_convolved = to_log_scale(psd_up_convolved)
psd_down_convolved = to_log_scale(psd_down_convolved)

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7))

visualize_min_distance = 10
visualize_max_distance = 20

x1 = visualize_min_distance * scene_scale
x2 = visualize_max_distance * scene_scale
x1p = x1 * B / (c * T) * 1e-6  # beat frequency in MHz
x2p = x2 * B / (c * T) * 1e-6  # beat frequency in MHz


# Up-Chirp
ax1.plot(distance, log_up_sample, label="Alg2. Single sample", color="C0", alpha=0.5)
ax1.plot(distance, log_up_avg, label="Alg2. Average", color="C0", linewidth=2)
ax1.plot(psd_distance, psd_up, label="Alg1. PSD", color="C1", linewidth=2, alpha=0.5)
ax1.plot(distance, psd_up_convolved, label="Alg1. PSD (conv)", color="C1", linewidth=2, linestyle="--")
ax1.set_title("Up-Chirp FFT Power (Log Scale)")
ax1.set_xlabel("Distance")
ax1.set_ylabel("Power (dB)")
ax1.set_xlim(x1, x2)
ax1.set_ylim(-30, 30)
ax1.grid(True)
ax1.legend()


# Down-Chirp
ax2.plot(distance, log_down_sample_flipped, label="Alg2. Single sample", color="C0", alpha=0.5)
ax2.plot(distance, log_down_avg_flipped, label="Alg2. Average", color="C0", linewidth=2)
ax2.plot(psd_distance, psd_down, label="Alg1. PSD", color="C1", linewidth=2, alpha=0.5)
ax2.plot(distance, psd_down_convolved, label="Alg1. PSD (conv)", color="C1", linewidth=2, linestyle="--")
ax2.set_title("Down-Chirp FFT Power (Log Scale)")
ax2.set_xlabel("Distance")
ax2.set_ylabel("Power (dB)")
ax2.set_xlim(x1, x2)
ax2.set_ylim(-30, 30)
ax2.grid(True)


# --- Secondary X axis ---
def idx_to_phys(x):
    return x1p + (x / (x2 - x1)) * (x2p - x1p)

def phys_to_idx(p):
    return (p - x1p) / (x2p - x1p) * (x2 - x1)

secax = ax1.secondary_xaxis('top', functions=(idx_to_phys, phys_to_idx))
secax.set_xlabel("Beat Frequency (MHz)")

secax = ax2.secondary_xaxis('top', functions=(idx_to_phys, phys_to_idx))
secax.set_xlabel("Beat Frequency (MHz)")

plt.tight_layout()

if not os.path.exists("../plots"):
    os.makedirs("../plots")
plt.savefig("../plots/%s_Alg1_Alg2_convergence.png" % scene_name)