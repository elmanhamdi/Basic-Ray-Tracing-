# CENG 488 Assignment8 by
# Elman Hamdi
# 240201036
# June 2021

from cameras import *
from utils import *
from objects import *
from scenes import *
from views import *
from lights import *
from qtui import *
import sys
import random
import time
import json
from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
import os

with open("file.json") as f:
    render_file = json.load(f)
    
render_settings = render_file['renderSettings']
WIDTH = render_settings['xres']
HEIGHT = render_settings['yres']
subsample = int(render_settings['subsample'])
nof_ambient_sample = int(render_settings['ambient_nof_sample'])
max_bounce = int(render_settings['max_bounce'])

# init scene
scene = Scene()
window = Window(HEIGHT, WIDTH)
camera = Camera(
    eye=Pos3d.list_to_pos(render_file['camera']['position']),
    center=Pos3d(0, 0, 0),
    window=window,
    window_distance=render_file['camera']['window_distance'],
)

for raw_sphere in render_file['sphere']:
    sphere = Sphere(
        radius=raw_sphere['radius'],
        center=Pos3d.list_to_pos(raw_sphere['position']),
        material=Material(
            color=Color(
                raw_sphere['color'][0],
                raw_sphere['color'][1],
                raw_sphere['color'][2],
            ),
            refractive=raw_sphere['refractive'],
            reflective=raw_sphere['reflective'],
        ),
    )
    scene.add(sphere)

for light in render_file['light']:
    l = Light(
    		light_type=light['type'], 
    		position=Pos3d.list_to_pos(light['position']),
              	color=Color(
              		light['color'][0], 
              		light['color'][1], 
              		light['color'][2]
              	), 
              	intensity=light['intensity'],
             )
    if l.light_type == 'dome':
        l.lightObj = Sphere(
        		radius=light['size'], 
        		center=l.position, 
        		material=Material(color=l.color),
        		)

    scene.add(l)

viewObj = View(	camera=camera, 
		scene=scene, 
		subsample=subsample, 
		ambient_sample=nof_ambient_sample,
               	max_bounce=max_bounce,
              )


if __name__ == "__main__":
    # setup a QApplication
    qApp = QApplication(sys.argv)
    qApp.setOrganizationName("CENG488 ")
    qApp.setOrganizationDomain("ELMAN")
    qApp.setApplicationName("Final Project")

    # setup main ui
    mainWindow = PyTraceMainWindow(qApp, WIDTH, HEIGHT)
    mainWindow.setupUi()
    mainWindow.show()

    # an example of writing to buffer
    mainTimer = QTimer()
    cpuCount = os.cpu_count()
    mainTimer.timeout.connect(mainWindow.timerBuffer(viewObj, num_thread=cpuCount - 1))
    mainTimer.start(2000)
    # enter event loop
    sys.exit(qApp.exec_())
