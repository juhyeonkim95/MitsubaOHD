import xml.etree.ElementTree as ET
import numpy as np
from utils.loader_utils import *


import copy
def set_node_key_value_boolean(node, key, value):
    key_node = node.find('*[@name="%s"]' % key)
    if key_node is None:
        child = ET.Element("boolean")
        child.set("name", key)
        child.set("value", "true" if value else "false")
        node.append(child)
    else:
        key_node.set("value", "true" if value else "false")


def set_node_key_value_float(node, key, value):
    key_node = node.find('*[@name="%s"]' % key)
    if key_node is None:
        child = ET.Element("float")
        child.set("name", key)
        child.set("value", str(value))
        node.append(child)
    else:
        key_node.set("value", str(value))

def create_element(tag, **kwargs):
    node = ET.Element(tag)
    for key, value in kwargs.items():
        if tag == "float" or tag == "integer":
            node.set(key, str(value))
        elif tag == "boolean" and key == "value":
            node.set(key, bool_to_str(value))
        else:
            node.set(key, value)
    return node


def add_scale_to_node(node, scale):
    if node.find("animation") != None:
        return
    if node.find("transform") == None:
        transform_node_temp = create_element("transform")
        transform_node_temp.set("name", "toWorld")
        node.append(transform_node_temp)
    
    transform_node = node.find("transform")
    scale_node = ET.Element("scale")
    scale_node.set("value", str(scale))
    transform_node.append(scale_node)
    


def add_animation_to_node(node, velocity, exposure_time):
    velocity_n = np.asarray(velocity, dtype=float) * exposure_time
    if node.find("animation") == None and node.find("transform") == None:
        transform_node_temp = create_element("transform")
        transform_node_temp.set("name", "toWorld")
        node.append(transform_node_temp)
    
    transform_node = node.find("transform")
    transform_node.attrib.pop("name")
    transform_node_new = copy.deepcopy(transform_node)
    translate_node = ET.Element("translate")
    translate_node.set("x", str(velocity_n[0]))
    translate_node.set("y", str(velocity_n[1]))
    translate_node.set("z", str(velocity_n[2]))
    transform_node_new.append(translate_node)

    transform_node.set("time", "0")
    transform_node_new.set("time", str(exposure_time))

    animation = ET.Element("animation")
    animation.set("name", "toWorld")
    animation.append(transform_node)
    animation.append(transform_node_new)

    node.remove(node.find("transform"))

    node.append(animation)
    