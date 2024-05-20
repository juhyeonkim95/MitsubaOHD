import numpy as np
from PIL import Image
import os
import matplotlib.pyplot as plt
import cv2


def rgb2luminance(img):
    return (0.2126 * img[:,:,0]) + (0.7152 * img[:,:,1]) + (0.0722 * img[:,:,2])

def ToneMap(c, limit):
    if len(c.shape) == 3 and c.shape[2] == 3:
        luminance = 0.3*c[:,:,0] + 0.6*c[:,:,1] + 0.1*c[:,:,2]
        luminance = np.dstack([luminance]*3)
    else:
        luminance = c
    col = c * 1.0 / (1.0 + luminance / limit)
    return col
def LinearToSrgb(c):
    kInvGamma = 1.0 / 2.2
    return np.power(c, kInvGamma)

def to_ldr_image(img):
    return LinearToSrgb(ToneMap(img, 1.5))


def save_rendered_img_ldr(img, folder, file_name):
    img_ldr = img ** (1.0 / 2.2)
    img_save = np.asarray(img_ldr * 255.0)
    img_save = np.clip(img_save, 0, 255)
    img_save = img_save.astype(np.uint8)
    im = Image.fromarray(img_save)
    if not os.path.exists(folder):
        os.makedirs(folder)
    im.save("%s/%s.png" % (folder, file_name))


def save_rendered_img_npy(img, folder, file_name):
    if not os.path.exists(folder):
        os.makedirs(folder)
    np.save("%s/%s.npy" % (folder, file_name), img)
    
def save_speed_image(image, output_path, filename, colorbar_also=False, velocity_range=5):
    cm = plt.get_cmap('RdBu')
    norm = plt.Normalize(-velocity_range, velocity_range)
    if len(image.shape) == 3:
        image = image[:, :, 0]

    image = cm(norm(image))
    image = image[:, :, 0:3]
    image = (image * 255.0).astype(np.uint8)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    cv2.imwrite(os.path.join(output_path, filename), image)


def save_tof_image(image, output_path, filename, 
    vmin_percentile=5, vmax_percentaile=95, 
    vmin=None, vmax=None, colorbar_also=False, crop_size=None):

    cm = plt.get_cmap('viridis')
    if vmin is None:
        vmin = np.percentile(image, vmin_percentile)
    if vmax is None:
        vmax = np.percentile(image, vmax_percentaile)
    norm = plt.Normalize(vmin, vmax)

    if len(image.shape) == 3:
        image = rgb2luminance(image)

    image = cm(norm(image))
    image = image[:, :, 0:3]
    image = (image * 255.0).astype(np.uint8)
    
    if crop_size is not None:
        start, end = crop_size
        image = image[start[1]:end[1], start[0]:end[0], :]

    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    cv2.imwrite(os.path.join(output_path, filename), image)


from skimage.util import random_noise
import skimage

def load_tof_image(file, add_noise=False):
    image = np.load(file)
    image = image - 1
    image = rgb2luminance(image)
    if add_noise:
        pass
        #noise = np.random.rand(*image.shape)
        #image += 1e-4 * noise
        # pass
        # noisy = random_noise(image, mode="poisson")
        # image = image * 0.95 + noisy * 0.05
        # image += np.random.poisson(image) * 1e-5
    return image

def load_velocity_image(file):
    return np.load(file)[:, :, 0] - 100.0
