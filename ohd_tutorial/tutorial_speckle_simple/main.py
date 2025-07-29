import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
import matplotlib as mat
mat.use("svg")
import os

def run(W):
    """ This function simulates simple speckle simulation for two-peaked scenario (Figure 7 in main paper). 
    It does three things
        (1) Evaluate normalized covariance matrix
        (2) Evaluate power spectrum  (show Alg1 and Alg2 converge to the same result)
        (3) Evaluate intensity histogram at specific frequency
    Args:
        W : width of the peak.
    """

    # Parameters
    N = 100  # Number of reflectors
    M = 1000  # Number of realizations (different configurations of random phases)
    N_T = 1024  # Number of time samples
    T = 1.0  # Total time duration
    sampling_rate = N_T / T  # Sampling rate in Hz

    f_min = 250  # Minimum frequency
    f_max = 350  # Maximum frequency
    f_1 = 275       # Peak 1 frequency
    f_2 = 325       # Peak 2 frequency

    signal_intensity = 1
    noise_intensity = 0.3

    # Time vector
    t = np.linspace(0, T, N_T, endpoint=False)

    # Random phases and frequencies for reflectors
    random_phases = np.random.uniform(-np.pi, np.pi, size=(M, N))  # Random phases 
    random_phases2 = np.random.uniform(-np.pi, np.pi, size=(M, N))  # Random phases 
    random_frequencies1 = np.random.uniform(f_1 - W, f_1 + W, size=(M, N))  # Random frequencies
    random_frequencies2 = np.random.uniform(f_2 - W, f_2 + W, size=(M, N))  # Random frequencies

    # Simulated signals
    intensities = []  # To store intensities at f_0 across all realizations
    
    for m in trange(M):
        # Create a signal for this realization
        signal = np.zeros(N_T, dtype=complex)
        signal_i = np.arange(N_T) / N_T
        for n in range(N):
            # Add contribution from each reflector with its random frequency and phase
            frequency1 = random_frequencies1[m, n]
            frequency2 = random_frequencies2[m, n]
            frequency = np.where(np.random.uniform() < 0.5, frequency1, frequency2)

            phase = random_phases[m, n]
            signal += (np.exp(1j * (2 * np.pi * frequency * t + phase)))

        # Normalize by sqrt(N)
        signal /= np.sqrt(N)
        # Perform FFT
        fft_result = np.fft.fft(signal)
        fft_result /= N_T
        fft_freqs = np.fft.fftfreq(N_T, d=1/sampling_rate)
        
        # add shot noise (shot noise also follows negative-exponential distribution)
        noise = np.random.exponential(scale=noise_intensity / N_T, size=N_T)
        intensities.append(np.abs(fft_result)**2 + noise)

    # Convert to numpy array for analysis
    intensities_all = np.array(intensities)
    total_power = np.sum(np.mean(intensities_all, axis=0))
    
    ################################
    # (1) Render Covariance Matrix #
    ################################
    if not os.path.exists("covariance_matrix"):
        os.makedirs("covariance_matrix")

    # Compute the normalized covariance matrix
    normalized_intensities_all = intensities_all / np.mean(intensities_all, axis=0, keepdims=True)
    cov_matrix = np.cov(normalized_intensities_all[:, f_min:f_max], rowvar=False)

    plt.figure(figsize=(4, 4))
    plt.imshow(cov_matrix, origin="lower", vmin=0.0, vmax=3.0, extent=[f_min, f_max, f_min, f_max])
    plt.xlabel("frequency")
    plt.ylabel("frequency")
    plt.title("Normalized Covariance Matrix")
    plt.tight_layout()
    plt.savefig("covariance_matrix/cov_matrix_%d.png" % W, dpi=600)
    plt.savefig("covariance_matrix/cov_matrix_%d.svg" % W, dpi=600)
    plt.close('all')

    #######################
    # (2) Render Spectrum #
    #######################
    if not os.path.exists("power_spectrum"):
        os.makedirs("power_spectrum")

    # (A) field sampling + random phase addition (Alg 2)
    intensities_all = intensities_all[:, f_min:f_max]
    intensity_mean = np.mean(intensities_all, axis=0)

    intensities_all_log = 10 * np.log10(intensities_all)
    intensity_mean_log = 10 * np.log10(intensity_mean)

    xf = fft_freqs[f_min:f_max]

    plt.figure(figsize=(4, 4))

    # plot two single measurements
    plt.plot(xf, intensities_all_log[1], color='C0', alpha=0.5, linewidth=1)
    plt.plot(xf, intensities_all_log[2], color='C1', alpha=0.5, linewidth=1)
    
    # plot ensemble averaged measurements
    plt.plot(xf, intensity_mean_log, color='red', linewidth=3, alpha=0.7)

    # (B) GT PSD (Alg 1 - no path sampling required for this simple scenario)
    # Define the x-axis range
    C = 10
    x = np.linspace(f_min, f_max, 100 * C)  # Adjust range as needed

    # Create the function
    # 1 if x is in [A, B], else 0
    y = np.where(((x >= f_1 - W) & (x <= f_1 + W)) | ((x >= f_2 - W) & (x <= f_2 + W)), 1.0, 0.0)
    y /= y.sum()
    y += noise_intensity / (N_T * C)
    y *= C # need to multiply because bin numbers are different!

    A = 1e-10
    ylog = 10 * np.log10(y + (A))
    plt.plot(x, ylog, color='black', linestyle='-', linewidth=1, alpha=0.5)

    # convolve with sinc square function.
    sinc_kernel = np.sinc((xf[:, None] - x[None, :]) / (xf[0] - xf[1])) ** 2
    sinc_kernel /= np.sum(sinc_kernel, axis=1, keepdims=True)

    ys = sinc_kernel @ y

    yslog = 10 * np.log10(ys + (A))
    plt.plot(xf, yslog, color='green', linestyle='-.', linewidth=1)


    plt.ylim(-45, -5)
    plt.xlabel("frequency")
    plt.ylabel("power (dB)")
    plt.title("Power Spectrum Density")
    plt.tight_layout()

    plt.savefig("power_spectrum/spectrum_%d.png" % W,dpi=600)
    plt.savefig("power_spectrum/spectrum_%d.svg" % W,dpi=600)
    plt.close('all')

    ########################
    # (3) Build Histogram  #
    ########################
    # only for W=9 case
    if not W == 9:
        return

    if not os.path.exists("intensity_histogram"):
        os.makedirs("intensity_histogram")

    target_freqs = [f_1-5, f_1, f_1+5, f_2-5, f_2, f_2+5]
    # target_freqs = [f_1, f_2]

    intensities_all = np.array(intensities)
    intensities_all /= np.sum(intensities_all, axis=0, keepdims=True)
    for target_freq in target_freqs:
        intensities = intensities_all[:, target_freq]

        # Plot histogram of intensities
        plt.figure(figsize=(8, 6))
        plt.hist(intensities, bins=50, density=True, alpha=0.7, color='blue', label='Simulated Intensity at Frequency %d' % target_freq)

        # Overlay the theoretical negative exponential distribution
        x = np.linspace(0, np.max(intensities), 1000)
        theoretical_pdf = (1 / np.mean(intensities)) * np.exp(-x / np.mean(intensities))
        plt.plot(x, theoretical_pdf, 'r-', lw=2, label='Negative Exponential (Theory)')

        # Add labels and legend
        plt.title('Intensity Distribution at Frequency %d' % target_freq, fontsize=14)
        plt.xlabel('Intensity (I)', fontsize=12)
        plt.ylabel('Probability Density', fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.savefig("intensity_histogram/freq_%d.png"%(target_freq))
        plt.close('all')

run(1)
run(3)
run(9)