import numpy as np

def get_velocity_dictionary():
    velocity_dictionary = {
        "default" : {
            "object": [0, 0, 0]
        },
        "doggler" : {
            "object": np.array([-2.5 * np.sin(np.rad2deg(20)), 2.5 * np.cos(np.rad2deg(20)), 0]) 
        },
        "cornell-box": {
            "ShortBox": [0, 0, 2],
            "TallBox": [0, 0, -2],
            "Floor": [0, 0.2, 0],
        },
        "living-room-2" : {
            "Mesh002": [0, 0, 2],
            "Mesh051": [0, 0, 2],
            "Mesh004": [0, 0, 2],
            "Mesh076": [0, 2, 0],
        },
        "road": {
            "car1":  [-2, 0, 0],
            "car2":  [2, 0, 0],
        }
    }
    return velocity_dictionary

def get_crop_pos(scene_name, crop_index=0):

    if "cornell-box" in scene_name:
        if "pos2" in scene_name:
            return (52, 122)
        elif "pos3" in scene_name:
            return (52, 64)
        return (32, 122)
        # return (31, 122)
    elif "living-room-2" in scene_name:
        # return (20, 33)
        return (12, 40)
    elif "road" in scene_name:
        return (110, 48)
    elif "rectangle" in scene_name:
        return (63, 65)
    elif "doggler" in scene_name:
        if crop_index == 0:
            return (245, 419)
        elif crop_index == 1:
            return (340, 419)
        elif crop_index == 2:
            return (378, 419)


def get_crop_dist_range(scene_name):
    if "cornell-box" in scene_name:
        return (11, 17)
    elif "living-room-2" in scene_name:
        return (2.5, 10)
    elif "road" in scene_name:
        return (5, 25)
    return (-1, -1)


def get_transient_min_max_range(scene_name):
    if "cornell-box" in scene_name:
        return (0, 40)
    elif "living-room-2" in scene_name:
        return (0, 20)
    elif "road" in scene_name:
        return (0, 100)
    elif "rectangle" in scene_name:
        return (0, 20)
    elif "doggler" in scene_name:
        return (0, 30)


def get_depth_min_max_range(scene_name):
    if "cornell-box" in scene_name:
        return (4, 10)
    elif "living-room-2" in scene_name:
        return (0, 4)
    elif "road" in scene_name:
        return (3, 20)
    # elif "rectangle" in scene_name:
    #     return (4, 12)
    
    
def get_vel_min_max_range(scene_name):
    if "cornell-box" in scene_name:
        return (-3, 3)
    elif "living-room-2" in scene_name:
        return (-3, 3)
    elif "road" in scene_name:
        return (-3, 3)
    
def get_default_fmcw_setting():
    configs = {}

    configs["radar"] = {
        'T': 7.33,
        'B': 0.15,
        'f_0': 77
    }

    configs["lidar_static"] = {
        'T': 1000,
        'B': 136,
        'wavelength': 1310
    }

    configs["lidar_dynamic"] = {
        'T': 10,
        'B': 1,
        'wavelength': 1550
    }

    return configs