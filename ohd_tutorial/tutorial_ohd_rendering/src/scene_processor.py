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

def add_fmcw_integrator(**kwargs):
    integrator_type = kwargs.get("integrator_type", "fmcwpsd")
    integrator = create_element("integrator", type=integrator_type)
    integrator.append(create_element("integer", name="maxDepth", value=kwargs.get("maxDepth", 4)))
    integrator.append(create_element("boolean", name="strictNormals", value="true"))
    integrator.append(create_element("integer", name="M", value=kwargs.get("M")))
    integrator.append(create_element("float", name="T", value=kwargs.get("T")))
    integrator.append(create_element("float", name="B", value=kwargs.get("B")))

    if "f_0" in kwargs:
        integrator.append(create_element("float", name="f_0", value=kwargs.get("f_0")))
    if "wavelength" in kwargs:
        integrator.append(create_element("float", name="wavelength", value=kwargs.get("wavelength")))

    if "psd" in integrator_type:    
        scene_scale = kwargs.get("scene_scale", 1)
        integrator.append(create_element("float", name="max_distance", value=kwargs.get("max_distance", 10) * scene_scale))
        integrator.append(create_element("float", name="min_distance", value=kwargs.get("min_distance", 0) * scene_scale))
        frames = kwargs.get("M") * 2
    else:
        frames = kwargs.get("M") * 4

    integrator.append(create_element("boolean", name="use_collimated", value=(kwargs.get("light_option", "point")=="collimated") ))
    
    scene_scale = kwargs.get("scene_scale", 1)

    if "fov_error" in kwargs:
        fov = kwargs.get("fov_error", 0.5)
        integrator.append(create_element("float", name="fov_error", value=fov))
    
    return integrator, frames


def process_scene_file_and_export(**kwargs):
    # load original scene file
    scene_files_path = kwargs.get("scene_files_path", "../scenes_original")
    scene_name = kwargs.get("scene_name", "cornell-box")
    scene_xml_name = kwargs.get("scene_xml_name", "scene.xml")
    scene_folder_path = os.path.join(scene_files_path, scene_name)
    scene_scale = kwargs.get("scene_scale", 1)

    # define output scene file
    new_scene_file_output_folder = kwargs.get("new_scene_file_output_folder", "../scenes_processed")
    output_folder = os.path.join(new_scene_file_output_folder, scene_name)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    xml_file_name = kwargs.get("output_xml_file_name", "scene.xml")
    
    # read scene file
    tree = ET.parse(os.path.join(scene_folder_path, scene_xml_name))
    root = tree.getroot()

    # prepare new scene file
    new_root = create_element("scene", version="0.6.0")

    # light configs
    light_option = kwargs.get("light_option", "point")

    # velocity dictionary
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

        # determine velocity
        velocity = velocity_dictionary_scene.get(ref_id, None)
        velocity = velocity_dictionary_scene.get(obj_file_name, velocity)
        if "road" in scene_name:
            if "001" in obj_file_name and "Material.001" not in obj_file_name:
                velocity = velocity_dictionary_scene.get("car1")
            elif "Plane" not in obj_file_name:
                velocity = velocity_dictionary_scene.get("car2") 

        # add velocity
        if velocity is not None and kwargs.get("add_velocity", True):
            velocity = np.asarray(velocity, dtype=float)
            if "fmcw" in integrator_type:
                add_velocity_to_node(node, velocity * scene_scale)
                # add_animation_to_node(node, velocity * 1e-6 * 0.2 * scene_scale, kwargs.get("T"))
            elif "velocity" in integrator_type:
                add_animation_to_node(node, velocity * scene_scale, 1e-6)
        
        # remove emitter?
        if light_option != "area" and node.find("emitter") is not None:
            if "cornell-box" in scene_name: 
                pass
            else:
                node.remove(node.find("emitter"))
                shapes.append(node)
        else:
            shapes.append(node)
        
    sensor_node = root.find("sensor")

    # scale scene
    if "road" not in scene_name and scene_scale != 1:
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
        if sensor_node.find("transform").find("matrix") != None:
            sensor_node.find("transform").find("matrix").set("value", new_value)
        else:
            sensor_node.remove(sensor_node.find('transform'))
            transform_new = create_element("transform", name="toWorld")
            matrix_new = create_element("matrix", value=new_value)
            transform_new.append(matrix_new)
            sensor_node.append(transform_new)

    if "road" in scene_name and scene_scale != 1:
        for t in ["x", "y", "z"]:
            v = sensor_node.find("transform").find("translate").get(t)
            sensor_node.find("transform").find("translate").set(t, str(float(v) * scene_scale))
    
    # (2) add point light if needed
    if light_option != "area":
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
    image_scale = kwargs.get("image_scale", 4)
    film_node.find('*[@name="width"]').set("value", str(width // image_scale))
    film_node.find('*[@name="height"]').set("value", str(height // image_scale))
    if film_node.find('*[@name="gamma"]') is not None:
        film_node.remove(film_node.find('*[@name="gamma"]'))
    if film_node.find('*[@name="banner"]') is not None:
        film_node.remove(film_node.find('*[@name="banner"]'))

    if kwargs.get("use_crop", False):
        crop_x, crop_y = get_crop_pos(scene_name, kwargs.get("crop_index", 0))
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
    integrator_type = kwargs.get("integrator_type", "fmcwpsd")
    if "fmcw" in integrator_type: 
        integrator, frames = add_fmcw_integrator(**kwargs)
    else:
        integrator, frames = add_depth_velocity_integrator(**kwargs)

    film_node.append(create_element("integer", name="frames", value=frames))
    new_root.append(integrator)

    # (FINAL) export scene
    new_tree = ET.ElementTree(new_root)
    ET.indent(new_tree, space="\t", level=0)
    new_tree.write(os.path.join(output_folder, xml_file_name))

    return os.path.join(output_folder, xml_file_name)