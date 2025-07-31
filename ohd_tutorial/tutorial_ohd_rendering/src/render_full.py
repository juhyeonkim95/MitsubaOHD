from scene_processor import *
import subprocess

def run_fmcwpsd_for_full_image(scene_name, scene_scale=10, max_depth=4):
    output_folder = os.path.join("../", "results")
    
    # set PSD min and max range (use distance unit instead of frequency for convenience)
    min_distance, max_distance = get_transient_min_max_range(scene_name)

    configs = {
        "scene_name": scene_name,
        "maxDepth": max_depth,
        "image_scale": 8,
        "add_velocity": True,
        "M": 1024,
        "B": 1,
        "T": 5,
        "wavelength": 1550,
        "scene_scale": scene_scale,
        "spp": 4096,
        "use_single_pixel": False,
        "integrator_type": "fmcwpsd",
        "light_option": "collimated", # point or collimated or area
        "min_distance": min_distance,
        "max_distance": max_distance,
    }

    output_file_name = "fmcwpsd"
    output_xml_file_name = output_file_name + ".xml"
    output_folder_name = os.path.join(output_folder, scene_name, "full")
    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)
    output_file_name = os.path.join(output_folder_name, output_file_name+".npy")

    new_file_path = process_scene_file_and_export(**configs, output_xml_file_name=output_xml_file_name)
    subprocess.run(["mitsuba", "-L", "error", "-o", output_file_name, new_file_path])


if __name__ == "__main__":
    run_fmcwpsd_for_full_images(scene_name="cornell-box-floor-specular", scene_scale=10, max_depth=4)