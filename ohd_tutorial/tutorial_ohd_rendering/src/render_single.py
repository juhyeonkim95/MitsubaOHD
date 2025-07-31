from scene_processor import *
import subprocess
from tqdm import trange

def run_for_single_pixel(scene_name, integrator, maxdepth=4, repeat=1, use_single_pixel=True):
    output_folder = os.path.join("../", "results")
    
    # set PSD min and max range (use distance unit instead of frequency for convenience)
    min_distance, max_distance = get_transient_min_max_range(scene_name)

    if integrator == "fmcwpsd":
        M = 4096
        spp= 4096 * 8
    else:
        M = 3072
        spp=256

    configs = {
        "scene_name": scene_name,
        "maxDepth": maxdepth,
        "image_scale": 8,
        "add_velocity": True,
        "M": M,
        "B": 1,
        "T": 5,
        "wavelength": 1550,
        "scene_scale": 10,
        "spp": spp,
        "use_single_pixel": use_single_pixel,
        "integrator_type": integrator,
        "light_option": "collimated", # point or collimated or area
        "min_distance": min_distance,
        "max_distance": max_distance,
    }

    output_folder_name = os.path.join(output_folder, scene_name, "single")
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)

    output_file_name = "%s_single_depth_%d" % (integrator, maxdepth)
    output_xml_file_name = output_file_name + ".xml"
    new_file_path = process_scene_file_and_export(**configs, output_xml_file_name=output_xml_file_name)

    # run renderer multiple times
    for i in trange(repeat):
        if repeat == 1:
            output_image_name = os.path.join(output_folder_name, output_file_name+".npy")
        else:
            if not os.path.exists(os.path.join(output_folder_name, output_file_name)):
                os.makedirs(os.path.join(output_folder_name, output_file_name))
            output_image_name = os.path.join(output_folder_name, output_file_name, "iter_%d.npy" % i)
        subprocess.run(["mitsuba", "-L", "error", "-o", output_image_name, new_file_path])


if __name__ == "__main__":
    use_single_pixel = True
    scene_name = "cornell-box-floor-specular"
    run_for_single_pixel(scene_name, "fmcwpsd", maxdepth=4, repeat=1, use_single_pixel=use_single_pixel)
    run_for_single_pixel(scene_name, "fmcwfield", maxdepth=4, repeat=1000, use_single_pixel=use_single_pixel)
    run_for_single_pixel(scene_name, "fmcwfield", maxdepth=2, repeat=1000, use_single_pixel=use_single_pixel)
    