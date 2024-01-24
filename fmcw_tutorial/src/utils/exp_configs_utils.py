def get_velocity_dictionary():
    velocity_dictionary = {
        "default" : {
            "object": [-10, 0, 0]
        },
        "cornell-box": {
            "ShortBox": [0, 0, 10],
            "TallBox": [0, 0, -10],
            "Floor": [0, 1, 0],
        },
        "bedroom" : {
            "Mesh065": [0, 0, 10],
            "Mesh067": [0, 0, 10],
            "Mesh033": [0, 0, 10],
            "Mesh025": [0, 0, 10],
            "Mesh052": [0, 0, 10]
        },
        "kitchen" : {
            "Mesh065": [0, 0, 10],
            "Mesh112": [0, 0, 10]
        },
        "veach-ajar" : {
            "Pot3": [0, 0, 10],
            "Pot2": [0, 0, 10],
            "Material": [0, 0, 10],
        },
        "living-room-2" : {
            "Mesh002": [0, 0, 2],
            "Mesh051": [0, 0, 2],
            "Mesh004": [0, 0, 2],
            "Mesh076": [0, 2, 0],
        },
        "soccer-ball" : {
            "mat-_Box_of_cerealBox_Box_of_cerealCOL_jpg": [-5, 0, 0],
            "mat-ball": [0, -10, 0],
            "Circle_002": [-5, 0, 0],
        },
        "road": {
            # "Carroceria_2_plano.003",
            # "Carroceria_2_plano",
            # "Logo_C_rculo_007",
            # "Luz_Tr_Plano_021",
            # "Marco_Parabrisad_Ad_Plano_006",
            # "Marco_Parabrisad_Ad_Plano_016",
            # "Marco_Ventana_2_Plano_012",
            # "Parabrisas_Tr_Plano_015",
            # "Pueratas_Plano_008",
            # "Rin_1_Ad_Plano_017",
            # "Rin_1_Tr_Plano_017",
            # "Ruedas_Ad_C_rculo_004",
            # "Ruedas_Tr_C_rculo_005",
        }
    }
    return velocity_dictionary

def get_crop_pos(scene_name):
    if "cornell-box" in scene_name:
        return (31, 122)
    elif "living-room-2" in scene_name:
        return (12, 40)
    elif "road" in scene_name:
        return (110, 48)
    elif "rectangle" in scene_name:
        return (63, 65)

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
        return (0, 40)
    elif "rectangle" in scene_name:
        return (0, 20)
    # elif scene_name == "bedroom":
    #     return (0, 20)
    # elif scene_name == "kitchen":
    #     return (0, 20)
    # elif scene_name == "cbox-bunny":
    #     return (0, 20)

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
        return (-1, 1)
    elif "road" in scene_name:
        return (-3, 3)
    
def get_default_fmcw_setting():
    configs = {}

    configs["radar"] = {
        'T': 7.33,
        'B': 0.15,
        'f_c': 77
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