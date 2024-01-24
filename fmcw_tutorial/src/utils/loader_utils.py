import re
import numpy as np

def str_to_bool(s):
    if s == "true" or s == "True":
        return True
    elif s == "false" or s == "False":
        return False
    return False

def str2floatarray(s):
    s = re.sub(' +', ' ', s)
    s = re.sub(', ', ',', s)
    fs = re.split(",| ", s)
    float_array = [float(x) for x in fs]
    return np.array(float_array, dtype=np.float64)


def str2_4by4mat(s):
    ms = str2floatarray(s)
    ms = ms.reshape((4, 4))
    return ms

def matrix44_to_str(matrix):
    flattened_array = np.array(matrix.transpose()).flatten()
    transform_str = ' '.join(map(str, flattened_array))
    return transform_str

def bool_to_str(value):
    return "true" if value else "false"

def float_array_to_str(float_array):
    s = ""
    for v in float_array:
        s+= "%.2f " % v
    return s