import os
from utils.xml_utils import *
from utils.exp_configs_utils import *
from utils.loader import *

def add_depth_velocity_integrator(**kwargs):
    integrator = create_element("integrator", type=kwargs.get("integrator_type"))
    if kwargs.get("integrator_type") == "depth":
        return integrator, 2
    else:
        return integrator, 1

def add_transient_integrator(**kwargs):
    integrator = create_element("integrator", type="transient")
    integrator.append(create_element("integer", name="maxDepth", value=kwargs.get("maxDepth", 4)))
    integrator.append(create_element("float", name="max_distance", value=kwargs.get("max_distance", 10)))
    integrator.append(create_element("float", name="min_distance", value=kwargs.get("min_distance", 0)))
    integrator.append(create_element("integer", name="bin", value=kwargs.get("bin", 1024)))
    integrator.append(create_element("boolean", name="strictNormals", value="true"))
    integrator.append(create_element("boolean", name="use_collimated", value=kwargs.get("use_collimated", False)))
    integrator.append(create_element("boolean", name="use_amplitude", value=kwargs.get("use_amplitude", False)))
    integrator.append(create_element("boolean", name="pdf_sqrt", value=kwargs.get("pdf_sqrt", False)))
    frames = kwargs.get("bin", 1024)

    scene_scale = kwargs.get("scene_scale", 1)
    if "fov_error" in kwargs:
        fov = kwargs.get("fov_error", 0.5) / scene_scale
        integrator.append(create_element("float", name="fov_error", value=fov))

        
    return integrator, frames

def add_linearity_integrator(**kwargs):
    output_mode = kwargs.get("output_mode", "single")
    if output_mode == "single":
        frames = 3
    elif output_mode == "transient":
        frames = 2 * (kwargs.get("M", 4) + 1)

    integrator = create_element("integrator", type="linearity")
    integrator.append(create_element("integer", name="maxDepth", value=kwargs.get("maxDepth", 4)))
    integrator.append(create_element("integer", name="target_depth", value=kwargs.get("maxDepth", 4)))
    integrator.append(create_element("boolean", name="strictNormals", value="true"))
    integrator.append(create_element("integer", name="M", value=kwargs.get("M", 4)))
    integrator.append(create_element("float", name="T", value=kwargs.get("T", 0.1)))
    integrator.append(create_element("string", name="output_mode", value=output_mode))
    integrator.append(create_element("string", name="spatial_correlation_mode", value=kwargs.get("spatial_correlation_mode", "ray_sampler")))
    return integrator, frames

def add_fmcw_integrator(**kwargs):
    integrator_type = kwargs.get("integrator_type", "fmcw")
    integrator = create_element("integrator", type=integrator_type)
    integrator.append(create_element("integer", name="maxDepth", value=kwargs.get("maxDepth", 4)))
    integrator.append(create_element("boolean", name="strictNormals", value="true"))
    integrator.append(create_element("integer", name="M", value=kwargs.get("M")))
    integrator.append(create_element("float", name="T", value=kwargs.get("T")))
    integrator.append(create_element("float", name="B", value=kwargs.get("B")))
    if "f_c" in kwargs:
        integrator.append(create_element("float", name="f_c", value=kwargs.get("f_c")))
    if "wavelength" in kwargs:
        integrator.append(create_element("float", name="wavelength", value=kwargs.get("wavelength")))
    if "invert" in kwargs:
        integrator.append(create_element("boolean", name="invert", value=kwargs.get("invert")))

    integrator.append(create_element("string", name="spatial_correlation_mode", value=kwargs.get("spatial_correlation_mode", "ray_sampler")))
    integrator.append(create_element("boolean", name="use_collimated", value=kwargs.get("use_collimated", False)))
    integrator.append(create_element("boolean", name="use_amplitude", value=kwargs.get("use_amplitude", False)))
    integrator.append(create_element("boolean", name="pdf_sqrt", value=kwargs.get("pdf_sqrt", False)))
    integrator.append(create_element("string", name="path_correlation", value=kwargs.get("path_correlation", "full")))
    
    scene_scale = kwargs.get("scene_scale", 1)
    if "fov_error" in kwargs:
        fov = kwargs.get("fov_error", 0.5) / scene_scale
        integrator.append(create_element("float", name="fov_error", value=fov))
    
    frames = kwargs.get("M")
    if integrator_type == "fmcwdiff" or integrator_type == "fmcwdynamic":
        frames *= 2
    
    return integrator, frames


def process_scene_file_and_export(**kwargs):
    # load original scene file
    scene_files_path = kwargs.get("scene_files_path", "../scenes_original")
    scene_name = kwargs.get("scene_name", "cornell-box")
    scene_xml_name = kwargs.get("scene_xml_name", "scene.xml")
    scene_folder_path = os.path.join(scene_files_path, scene_name)
    scene_scale = kwargs.get("scene_scale", 1)
    

    # define output scene file
    output_scene_files_path = kwargs.get("output_scene_files_path", "../scenes_processed")
    if scene_scale != 1:
        output_folder = os.path.join(output_scene_files_path, scene_name+"_scale_%d" % scene_scale)
    else:
        output_folder = os.path.join(output_scene_files_path, scene_name)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    xml_file_name = kwargs.get("output_xml_file_name", "scene.xml")
    
    # read scene file
    tree = ET.parse(os.path.join(scene_folder_path, scene_xml_name))
    root = tree.getroot()

    # prepare new scene file
    new_root = create_element("scene", version="0.6.0")

    # other configs
    force_collocated_point_light = kwargs.get("force_collocated_point_light", True)
    velocity_dictionary = kwargs.get("velocity_dictionary", get_velocity_dictionary())
    if "cornell-box" in scene_name:
        velocity_dictionary_scene = velocity_dictionary.get("cornell-box")
    else:
        velocity_dictionary_scene = velocity_dictionary.get(scene_name, velocity_dictionary.get("default"))

    # (1) find shape and add velocity if needed
    shapes = []
    integrator_type = kwargs.get("integrator_type", "linearity")
    for node in root.findall('shape'):
        add_scale_to_node(node, scene_scale)

        # reference for velocity dict
        ref_id = None
        if node.find("ref") is not None:
            ref_id = node.find("ref").attrib["id"]
        obj_file_name = None
        if node.attrib["type"] == "obj" or node.attrib["type"] == "ply":
            obj_file_name = node.find('*[@name="filename"]').attrib["value"]
            obj_file_name = os.path.splitext(os.path.basename(obj_file_name))[0]

        velocity = velocity_dictionary_scene.get(ref_id, None)
        velocity = velocity_dictionary_scene.get(obj_file_name, velocity)
        if "road" in scene_name:
            if "001" in obj_file_name and "Material.001" not in obj_file_name:
                velocity = [-10, 0, 0]
            elif "plane" not in obj_file_name:
                velocity = [10, 0, 0]

        if velocity is not None and kwargs.get("add_velocity", True):
            velocity = np.asarray(velocity, dtype=float)
            if "linear" in integrator_type:
                add_animation_to_node(node, velocity * 0.2 * scene_scale, kwargs.get("T"))
                # add_animation_to_node(node, velocity, 1.0)
            elif "fmcw" in integrator_type:
                add_animation_to_node(node, velocity * 1e-6 * 0.2 * scene_scale, kwargs.get("T"))
            elif "velocity" in integrator_type:
                add_animation_to_node(node, velocity * 0.2 * scene_scale, 1e-6)
        
        if force_collocated_point_light and node.find("emitter") is not None:
            if "cornell-box" in scene_name: 
                pass
            else:
                node.remove(node.find("emitter"))
                shapes.append(node)
        else:
            shapes.append(node)
        
    sensor_node = root.find("sensor")

    if "road" not in scene_name and scene_scale != 1:
        # transform_node = sensor_node.find("transform")
        transform_matrix = load_value(sensor_node, "toWorld")
        transform_matrix[3][0] *= scene_scale
        transform_matrix[3][1] *= scene_scale
        transform_matrix[3][2] *= scene_scale
        temp = np.asarray(transform_matrix).T.flatten()
        new_value = ""
        for i, t in enumerate(temp):
            new_value += str(t)
            if i < len(temp) - 1:
                new_value += " "
        sensor_node.find("transform").find("matrix").set("value", new_value)

    if "road" in scene_name and scene_scale != 1:
        for t in ["x", "y", "z"]:
            v = sensor_node.find("transform").find("translate").get(t)
            sensor_node.find("transform").find("translate").set(t, str(float(v) * scene_scale))

    # (2) add point light if needed
    if force_collocated_point_light:
        point_light = create_element("emitter", type="point")
        point_light.append(sensor_node.find("transform"))
        point_light.append(create_element("rgb", name="intensity", value="100"))
        new_root.append(point_light)

    # (3) modify file paths
    for bsdf in root.findall('bsdf'):
        for filename in bsdf.findall('.//*[@name="filename"]'):
            filename.set("value", "%s/%s" % (scene_folder_path, filename.get("value")) )
        new_root.append(bsdf)
    for shape in shapes:
        for filename in shape.findall('.//*[@name="filename"]'):
            filename.set("value", "%s/%s" % (scene_folder_path, filename.get("value")))
        new_root.append(shape)

    # (4) add sensor node
    # (4-1) file node
    film_node = sensor_node.find("film")
    film_node.set("type", "mfilm")
    film_node.find('*[@name="fileFormat"]').set("value", "numpy")
    film_node.find('*[@name="pixelFormat"]').set("value", "luminance")
    width = load_value(film_node, "width")
    height = load_value(film_node, "height")
    scale = kwargs.get("scale", 4)
    film_node.find('*[@name="width"]').set("value", str(width // scale))
    film_node.find('*[@name="height"]').set("value", str(height // scale))
    if film_node.find('*[@name="gamma"]') is not None:
        film_node.remove(film_node.find('*[@name="gamma"]'))
    if film_node.find('*[@name="banner"]') is not None:
        film_node.remove(film_node.find('*[@name="banner"]'))

    if kwargs.get("use_crop", False):
        crop_x, crop_y = get_crop_pos(scene_name)
        film_node.append(create_element("integer", name="cropOffsetX", value=crop_x))
        film_node.append(create_element("integer", name="cropOffsetY", value=crop_y))
        film_node.append(create_element("integer", name="cropWidth", value=1))
        film_node.append(create_element("integer", name="cropHeight", value=1))
    
    
    # (4-2) sampler
    sampler_node = sensor_node.find("sampler")
    sampler_node.set("type", "independent")
    sampler_node.find('*[@name="sampleCount"]').set("value", (str(kwargs.get("spp", 512))))
    new_root.append(sensor_node)

    # (5) add integrator
    if integrator_type == "linearity":
        integrator, frames = add_linearity_integrator(**kwargs)
    elif integrator_type == "transient":
        integrator, frames = add_transient_integrator(**kwargs)
    elif integrator_type == "depth" or integrator_type == "velocity":
        integrator, frames = add_depth_velocity_integrator(**kwargs)
    else:
        integrator, frames = add_fmcw_integrator(**kwargs)
    film_node.append(create_element("integer", name="frames", value=frames))
    new_root.append(integrator)


    # (FINAL) export scene
    new_tree = ET.ElementTree(new_root)
    ET.indent(new_tree, space="\t", level=0)
    new_tree.write(os.path.join(output_folder, xml_file_name))

    return os.path.join(output_folder, xml_file_name)



