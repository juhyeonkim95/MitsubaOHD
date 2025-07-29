import numpy as np
import math


def lerp(a, b, f):
    return a * (1-f) + b * f

def normalize(mat):
    return mat / np.linalg.norm(mat)


def length2(vec):
    return vec.dot(vec)


def length(vec):
    return math.sqrt(vec.dot(vec))
