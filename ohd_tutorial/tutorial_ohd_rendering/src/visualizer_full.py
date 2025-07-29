import numpy as np
import matplotlib.pyplot as plt
from utils.exp_configs_utils import *

scene_name = "cornell-box-floor-specular"
min_distance, max_distance = get_transient_min_max_range(scene_name)
scene_scale = 10

# --- Load and Flip Data ---
data = np.load("../results/%s/full/fmcwpsd.npy" % scene_name)
data = np.flipud(data)

H, W, double_M = data.shape
M = double_M // 2

# FMCW constants
T = 5e-6
B = 1e9
c = 3e8

# --- x-axis settings ---
x1 = min_distance * scene_scale
x2 = max_distance * scene_scale
x1p = x1 * B / (c * T) * 1e-6  # beat frequency in MHz
x2p = x2 * B / (c * T) * 1e-6  # beat frequency in MHz
x_range = np.linspace(x1, x2, M)

eps = 1e-6
log_scale = True
current_coords = [0, 0]

# --- Precompute average image ---
avg_image = data.mean(axis=2)

# --- Setup plot ---
fig, (ax_img, ax_plot) = plt.subplots(1, 2, figsize=(12, 6))

# Left heatmap
im = ax_img.imshow(avg_image, cmap='viridis', origin='lower')
ax_img.set_title("Average of PSD (no meaning. just for visualization)")
cursor_dot, = ax_img.plot([], [], 'r+', ms=8)

# Right spectrum plot
plot1, = ax_plot.plot([], [], label='Up-Chirp', lw=2)
plot2, = ax_plot.plot([], [], linestyle='dotted', label='Down-Chirp', lw=2)

ax_plot.set_xlim(x1, x2)
ax_plot.set_xlabel("Distance (m)")
ax_plot.set_ylabel("Power (dB)" if log_scale else "Power")
ax_plot.set_title("Selected Pixel Spectrum")
ax_plot.set_yscale('linear')
ax_plot.legend()

# --- Toggle display text ---
toggle_text = ax_plot.text(
    0.98, 0.85,
    "Power Scale: dB (press 'j' to toggle)",
    transform=ax_plot.transAxes,
    fontsize=10,
    ha='right',
    va='top',
    bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round')
)

# --- Secondary X axis ---
def idx_to_phys(x):
    return x1p + (x / (x2 - x1)) * (x2p - x1p)

def phys_to_idx(p):
    return (p - x1p) / (x2p - x1p) * (x2 - x1)

secax = ax_plot.secondary_xaxis('top', functions=(idx_to_phys, phys_to_idx))
secax.set_xlabel("Beat Frequency (MHz)")

# --- Update Plot ---
def update_plot(x, y):
    global current_coords
    current_coords = [x, y]

    if not (0 <= x < W and 0 <= y < H):
        return

    vec = data[y, x]
    y1 = vec[:M]
    y2 = vec[M:]

    if log_scale:
        y1p = 10 * np.log10(np.clip(y1, eps, None))
        y2p = 10 * np.log10(np.clip(y2, eps, None))
    else:
        y1p = y1
        y2p = y2

    if np.all(np.isneginf(y1p)) and np.all(np.isneginf(y2p)):
        plot1.set_data([], [])
        plot2.set_data([], [])
        ax_plot.set_ylim(-100, 0)
        ax_plot.set_title(f"No valid dB values at ({x}, {y})")
    else:
        plot1.set_data(x_range, y1p)
        plot2.set_data(x_range, y2p)
        ymin = min(np.min(y1p), np.min(y2p))
        ymax = max(np.max(y1p), np.max(y2p))
        ax_plot.set_ylim(ymin, ymax)
        ax_plot.set_title(f"Selected Pixel Spectrum at ({x}, {y})")

    fig.canvas.draw_idle()


# --- Handlers ---
def on_mouse_move(event):
    if event.inaxes != ax_img:
        return
    x = int(event.xdata + 0.5)
    y = int(event.ydata + 0.5)
    cursor_dot.set_data([x], [y])
    update_plot(x, y)

def on_key_press(event):
    global log_scale
    if event.key == 'j':
        log_scale = not log_scale
        ax_plot.set_ylabel("Power (dB)" if log_scale else "Power")
        toggle_text.set_text("Power Scale: dB (press 'j' to toggle)" if log_scale
                             else "Power Scale: Linear (press 'j' to toggle)")
        update_plot(*current_coords)


# --- Connect ---
fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
fig.canvas.mpl_connect('key_press_event', on_key_press)

plt.tight_layout()
plt.show()
