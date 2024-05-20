import os
import numpy as np
import matplotlib.pyplot as plt
from utils.image_utils import *
from utils.common_path import *
import itertools
from scipy.fft import fft, ifft, fftfreq, fftshift
from utils.exp_configs_utils import *



def export_depth_map(**kwargs):
    # all of the different options
    scene_names = kwargs.get("scene_names")
    maxDepths = kwargs.get("maxDepths")
    light_option_names = kwargs.get("light_option_names")
    scales = kwargs.get("scales")
    data_options = kwargs.get("data_options")
    cropOffsetX = kwargs.get("cropOffsetX", -1)
    cropOffsetY = kwargs.get("cropOffsetY", -1)
    cropWidth = kwargs.get("cropWidth", -1)
    cropHeight = kwargs.get("cropHeight", -1)
    spps = kwargs.get("spps")
    file_names = kwargs.get("file_names")
    scene_scale = kwargs.get("scene_scale",1)

    # FMCW related
    T = kwargs.get("T", 10) * 1e-6
    B = kwargs.get("B", 1) * 1e9
    wavelength = kwargs.get("wavelength", 1550) * 1e-9
    c  = 3e8
    f_c = c / wavelength

    image_output_folder = kwargs.get("image_output_folder")
    
    images = {}
    first_image = None
    bin_size = 2048
    M = 4096

    for scene_name in scene_names:
        scene_name_org = scene_name
        if scene_scale != 1:
            scene_name = "%s_scale_%d" % (scene_name, scene_scale)

        fp1 = None

        for i in range(len(file_names)):
            output_xml_file_name = file_names[i]
            exp_output_folder = os.path.join(image_output_folder, data_options[i], scene_name)
            output_file_name = os.path.join(exp_output_folder, output_xml_file_name)
            image = np.load(output_file_name+".npy")
            
            perform_depthmap_calc = False

            if data_options[i] == "fmcw_dynamic":
                channels = image.shape[2]
                yf1 = fft(image[:,:,0:channels // 2], axis=-1)
                magnitude1 = np.abs(yf1)
                magnitude_max_idx1 = np.argmax(magnitude1, axis=-1)

                yf2 = fft(image[:,:,channels//2:channels], axis=-1)
                magnitude2 = np.abs(yf2)
                magnitude_max_idx2 = np.argmax(magnitude2, axis=-1)
                
                xf = fftfreq(M, T) * M
                fp1 = np.abs(xf[magnitude_max_idx1])
                fp2 = np.abs(xf[magnitude_max_idx2])
                perform_depthmap_calc = True
            else:
                yf = fft(image, axis=-1)
                magnitude = np.abs(yf)
                magnitude_max_idx = np.argmax(magnitude, axis=-1)
                xf = fftfreq(M, T) * M
                fp = np.abs(xf[magnitude_max_idx])

                if (i % 2) == 0:
                    fp1 = fp
                    perform_depthmap_calc = False
                else:
                    fp2 = fp
                    perform_depthmap_calc = True
                
            if perform_depthmap_calc:
                f_R = ((fp1 + fp2) * 0.5)
                depth_map = f_R * c * T / B * 0.5

                f_D = ((fp1 - fp2) * 0.5)
                velocity_map = f_D * c / f_c * 0.5
                
                depthimage_output_folder = os.path.join(kwargs.get("depthmap_output_folder"), scene_name)
                if not os.path.exists(depthimage_output_folder):
                    os.makedirs(depthimage_output_folder)

                dmin, dmax = get_depth_min_max_range(scene_name_org)
                plt.imshow(depth_map, vmin=dmin * scene_scale, vmax=dmax * scene_scale)
                plt.axis('off')
                plt.savefig(os.path.join(depthimage_output_folder, "depth_%s.png" % output_xml_file_name), bbox_inches='tight', pad_inches=0, dpi=600)
                plt.close('all')

                fig = plt.figure()
                vmin, vmax = get_vel_min_max_range(scene_name_org)
                plt.imshow(velocity_map, vmin=vmin * scene_scale, vmax=vmax * scene_scale, cmap="RdBu")
                plt.axis('off')
                plt.savefig(os.path.join(depthimage_output_folder, "vel_%s.png" % output_xml_file_name), bbox_inches='tight', pad_inches=0, dpi=600)
                plt.close('all')

        
if __name__ == "__main__":
    scene_names = [
        "cornell-box",
    ]
    data_options = ["fmcw_dynamic"]
    
    file_names = [
        "fmcw_dynamic_duration_10_depth_2_scale_8_spp_1_collimated_amplitude",
    ]
    
    configs = {
        "scene_names": scene_names,
        "file_names": file_names,
        "image_output_folder": os.path.join(PROJECT_PATH, "results"),
        "depthmap_output_folder" : "../results_fft/depth_velocity_map_recon",
        "scene_scale":10,
        "data_options": data_options
    }
    export_depth_map(
        **configs
    )