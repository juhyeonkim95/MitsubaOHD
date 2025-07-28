import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
mat.use('svg')
import torch
import torch.nn.functional as F

def generate_gaussian_surface(sigma_H, l, W, N, seed):
    """
    Generate a 2D Gaussian random surface.

    Parameters:
        sigma_H (float): RMS height of the surface.
        l (float): Correlation length.
        W (float): Patch size (physical size of the surface).
        N (int): Number of points along one dimension (creates an NxN array).

    Returns:
        np.ndarray: NxN array containing height data.
    """
    # Define grid spacing and frequency domain grid
    dx = W / N
    kx = np.fft.fftfreq(N, d=dx)
    ky = np.fft.fftfreq(N, d=dx)
    kx, ky = np.meshgrid(kx, ky)
    k = np.sqrt(kx**2 + ky**2)

    # Define power spectral density (PSD) for Gaussian surface
    PSD = np.exp(-k**2 * l**2 / 2)
    np.random.seed(seed)

    # Generate random phase and amplitude
    random_phase = np.exp(1j * 2 * np.pi * np.random.rand(N, N))
    amplitude = np.sqrt(PSD)

    # Create the Fourier transform of the surface
    surface_ft = amplitude * random_phase

    # Transform to spatial domain and scale by RMS height
    surface = np.fft.ifft2(surface_ft).real

    # Normalize RMS height to sigma_H
    surface *= sigma_H / np.sqrt(np.mean(surface**2))

    return surface


def sample_heightmap_cuda(surface, W, A, x, y, M, rx, ry):
    """
    Extract a sampled heightmap from the generated surface using CUDA grid_sample.

    Parameters:
        surface (np.ndarray): The generated heightmap (NxN array).
        W (float): Patch size of the full surface (e.g., in micrometers).
        A (float): Size of the area to sample (e.g., in micrometers).
        x (float): X-coordinate of the center of the sampling area (e.g., in micrometers).
        y (float): Y-coordinate of the center of the sampling area (e.g., in micrometers).
        M (int): Size of the output sampled heightmap (MxM array).

    Returns:
        np.ndarray: Interpolated MxM sampled heightmap.
    """
    N = surface.shape[0]  # Grid size of the full surface

    # Normalize coordinates for grid_sample (range [-1, 1])
    normalized_x = (x + W / 2) / W * 2 - 1
    normalized_y = (y + W / 2) / W * 2 - 1
    normalized_A = A / W * 2

    # Generate normalized sampling grid
    grid_x = torch.linspace(normalized_x - normalized_A / 2, normalized_x + normalized_A / 2, M, device='cuda').double()
    grid_y = torch.linspace(normalized_y - normalized_A / 2, normalized_y + normalized_A / 2, M, device='cuda').double()
    grid_y, grid_x = torch.meshgrid(grid_y, grid_x, indexing='ij')
    # grid = torch.stack((grid_x, grid_y), dim=-1).unsqueeze(0)  # Shape: (1, M, M, 2)

    perturbation_scale = 1.0
    # Add random perturbations to the sampling grid
    # perturbation_x = rx  # Range [-1, 1]
    # perturbation_y = ry  # Range [-1, 1]
    # perturbation_x *= (perturbation_scale * 2 / M)  # Scale based on grid spacing
    # perturbation_y *= (perturbation_scale * 2 / M)  # Scale based on grid spacing
    perturbed_grid_x = grid_x + rx * (2 / M)
    perturbed_grid_y = grid_y + ry * (2 / M)

    grid = torch.stack((perturbed_grid_x, perturbed_grid_y), dim=-1).unsqueeze(0)  # Shape: (1, M, M, 2)


    # Use grid_sample for interpolation
    sampled_tensor = F.grid_sample(surface, grid, mode='bilinear')

    # Convert back to numpy
    sampled_surface = sampled_tensor.squeeze()

    # Local coordinates in the sampling area
    local_x = torch.linspace(-A / 2, A / 2, M, device='cuda')
    local_y = torch.linspace(-A / 2, A / 2, M, device='cuda')
    local_y, local_x = torch.meshgrid(local_y, local_x, indexing='ij')

    return torch.stack([local_x, local_y, sampled_surface], dim=-1)

def plot_surface(surface, filename, W=400, **kwargs):
    # Plot the surface
    plt.figure(figsize=(8, 6))
    plt.imshow(surface, origin='lower', cmap='terrain', extent=[-W/2, W/2,-W/2,W/2], **kwargs)
    plt.colorbar(label='Height')
    plt.title("2D Gaussian Random Surface")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(filename)
    plt.close('all')
