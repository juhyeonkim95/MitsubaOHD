import numpy as np
import matplotlib.pyplot as plt
from utils.exp_configs_utils import *
import matplotlib
matplotlib.use("svg")
import os

# --- Configuration ---
scene_name = "cornell-box-floor-specular"
min_distance, max_distance = get_transient_min_max_range(scene_name)
scene_scale = 10

# --- Load Data ---
data = np.load(f"../results/{scene_name}/full/fmcwpsd.npy")

H, W, double_M = data.shape
M = double_M // 2

data_up = data[:, :, 0:M]
data_down = data[:, :, M:2*M]

# --- Constants ---
T = 5e-6
B = 1e9
c = 3e8
f_c = c / 1550e-9  # center frequency for 1550 nm

# --- x-axis for distance ---
x1 = min_distance * scene_scale
x2 = max_distance * scene_scale
xs = np.linspace(x1, x2, M)

# --- Extract peak positions ---
power_max_idx_up = np.argmax(data_up, axis=-1)
power_max_up = xs[power_max_idx_up]

power_max_idx_down = np.argmax(data_down, axis=-1)
power_max_down = xs[power_max_idx_down]

# --- Compute depth and velocity ---
distance_map = (power_max_up + power_max_down) * 0.5 * 0.5
velocity_map = (power_max_down - power_max_up) * 0.5 * 0.5
velocity_map /= (T * f_c / B)

# --- Plot and save ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

dmin, dmax = get_depth_min_max_range(scene_name)
vmin, vmax = get_vel_min_max_range(scene_name)

im1 = ax1.imshow(distance_map, cmap="viridis", vmax= dmax* scene_scale, vmin=dmin * scene_scale)
ax1.set_title("Distance (m)")
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

im2 = ax2.imshow(velocity_map, cmap="RdBu", vmax= vmax* scene_scale, vmin=vmin * scene_scale)
ax2.set_title("Velocity (m/s)")
plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

plt.tight_layout()

if not os.path.exists("../plots"):
    os.makedirs("../plots")
plt.savefig("../plots/%s_distance_velocity_reconstruction.png" % scene_name)