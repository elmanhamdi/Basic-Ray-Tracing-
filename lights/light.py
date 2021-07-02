# CENG 488 Assignment8 by
# Elman Hamdi
# 240201036
# June 2021

from utils import *
from objects import  *

#type can be 'dome' or 'point'
class Light:
    def __init__(self, light_type, position=Pos3d(0, 0, 0), color=[1, 1, 1, 1], intensity=1, size = None, lightObj = None):
        self.light_type = light_type
        self.position = position
        self.color = color
        self.intensity = intensity
        self.size = None
        self.lightObj = lightObj


class Color:
    def __init__(self, r, g, b, a):
        self.r = r
        self.g = g
        self.b = b
        self.a = a

    def __str__(self):
        return '[' + self.r + ', ' + self.g + ', ' + self.b + ', ' + self.a + ']'
