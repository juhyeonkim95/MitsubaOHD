import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from tqdm import tqdm, trange
import numpy as np
from scipy.fft import fft, ifft, fftfreq, fftshift
import matplotlib
matplotlib.use('Svg')
import tamaas as tm


def save_heightmap(heightmap, filename, wave_length, actual_patch_size):
    unit_length = wave_length
    
    # Calculate the number of units on each axis
    num_units_x = int(actual_patch_size / unit_length)
    num_units_y = int(actual_patch_size / unit_length)

    # Generate tick positions for the desired unit length scaling
    x_ticks = np.arange(-num_units_x // 2, num_units_x // 2 + 1) * unit_length
    y_ticks = np.arange(-num_units_y // 2, num_units_y // 2 + 1) * unit_length

    # Filter tick labels to show only multiples of 10
    x_labels = [int(x / unit_length) if x % (10 * unit_length) == 0 else '' for x in x_ticks]
    y_labels = [int(y / unit_length) if y % (10 * unit_length) == 0 else '' for y in y_ticks]

    # Plot the image
    plt.imshow(heightmap / wave_length, cmap='terrain', extent=[-actual_patch_size / 2, actual_patch_size / 2,
                                -actual_patch_size / 2, actual_patch_size / 2])

    # Set the ticks and labels
    plt.xticks(ticks=x_ticks, labels=x_labels)
    plt.yticks(ticks=y_ticks, labels=y_labels)

      # 'terrain' colormap is suitable for heightmaps
    plt.colorbar()  # Add a colorbar to show height values
    plt.title('Heightmap')
    plt.savefig("%s.png"% filename)
    plt.close('all')

    # acf = tm.Statistics2D.computeAutocorrelation(heightmap)
    # acf = np.fft.fftshift(acf)

    # plt.figure()
    # plt.imshow(acf)
    # plt.savefig("%s_acf.png"% filename)
    # plt.close('all')
