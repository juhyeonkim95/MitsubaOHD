from utils.loader_simple import *


def load_value(node, name, default=None, key="value"):
    """
    Load value from xml node by name
    :param node: xml node
    :param name: name of value
    :param default: default value
    :param key: key of load target. default is "value"
    :return: loaded value
    """
    child_node = node.find('*[@name="%s"]' % name)
    if child_node is not None:
        tag = child_node.tag
        # scalar / vector / matrix
        if tag == "integer":
            return load_int(child_node, key)
        if tag == "float":
            return load_float(child_node, key)
        elif tag == "boolean":
            return load_boolean(child_node, key)
        elif tag == "string":
            return load_string(child_node, key)
        elif tag == "point" or tag == "vector":
            return load_vector(child_node, key, default=0)
        elif tag == "transform":
            return load_transform(child_node)
    else:
        return default