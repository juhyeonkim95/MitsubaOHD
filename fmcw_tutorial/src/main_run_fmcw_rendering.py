from program_runner import *
import itertools
import subprocess
import matplotlib.pyplot as plt
from utils.common_path import *
import configargparse

def run(**kwargs):
    scene_names = kwargs.get("scene_names")
    maxDepths = kwargs.get("maxDepths")
    output_folder = kwargs.get("output_folder")

    if kwargs.get("use_collimated", False):
        light_option_name = "collimated"
    elif kwargs.get("force_collocated_point_light", True):
        light_option_name = "point"
    else:
        light_option_name = "area"
    if kwargs.get("use_amplitude", False):
        light_option_name += "_amplitude"
    if kwargs.get("pdf_sqrt", False):
        light_option_name += "_pdfsqrt"

    scales = kwargs.get("scales")
    T=kwargs.get("T", 10)
    scene_scale=kwargs.get("scene_scale", 1)
    spp=kwargs.get("spp")
    path_correlation = kwargs.get("path_correlation")
    fov = kwargs.get("fov_error")
    configs = list(itertools.product(scene_names, maxDepths, scales))
    
    for config in configs:
        scene_name, maxDepth, scale = config

        def run(output_xml_file_name, invert=False, has_fov=False):
            kwargs2 = {
                "scene_name": scene_name,
                "output_xml_file_name": output_xml_file_name + ".xml",
                "maxDepth": maxDepth,
                "scale": scale,
                "invert": invert,
                **kwargs
            }

            new_file_path = process_scene_file_and_export(**kwargs2)
            exp_output_folder = os.path.join(output_folder, scene_name)
            if scene_scale != 1:
                exp_output_folder = os.path.join(output_folder, scene_name+"_scale_%d" % scene_scale)

            if not os.path.exists(exp_output_folder):
                os.makedirs(exp_output_folder)
            output_file_name = os.path.join(exp_output_folder, output_xml_file_name)
            if os.path.isfile(output_file_name+".npy"):
                print(output_file_name, "EXIST")
            else:
                print(output_file_name, "RENDER")
            subprocess.run(["mitsuba", "-L", "error", "-o", output_file_name, new_file_path])

        if kwargs.get("runtype") == "run_dynamic_single" or kwargs.get("runtype") == "run_dynamic":
            output_xml_file_name_temp = "fmcw_dynamic_duration_%d_depth_%d_scale_%d_spp_%d_%s" % ((int)(T), maxDepth, scale, spp, light_option_name)
            run(output_xml_file_name_temp)
        elif kwargs.get("runtype") == "run_static_single" or kwargs.get("runtype") == "run_static":
            output_xml_file_name_temp = "fmcw_static_duration_%d_depth_%d_scale_%d_spp_%d_%s" % ((int)(T), maxDepth, scale, spp, light_option_name)
            run(output_xml_file_name_temp)

def run_dynamic():
    scene_names = [
        "cornell-box",
    ]
    maxDepths = [2]
    scales = [8]

    output_folder = os.path.join(PROJECT_PATH, "results/fmcw_dynamic")

    configs = {
        "scene_names": scene_names,
        "maxDepths": maxDepths,
        "scales": scales,
        "output_folder": output_folder,
        "spp": 1,
        "add_velocity": True,
        "M": 4096,
        "B": 1,
        "T": 10,
        "wavelength": 1550,
        "scene_scale": 10,
        "integrator_type": "fmcwdynamic",
        "runtype": "run_dynamic",
        "use_amplitude":True,
        "pdf_sqrt":False,
    }

    run(**configs, use_collimated=True)


if __name__ == "__main__":
    run_dynamic()