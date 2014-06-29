# A short script showing how to use the mesh_kit_2D to create a mesh of quarter_annulus,
# and plots the resulting mesh using matplotlib.
#
# Created 29.6 2014

import numpy as np
import matplotlib.pyplot as plt
import mesh_kit_2D

mesh = mesh_kit_2D.quarter_annulus_2D(0.01)
points = np.array(mesh.points)
tri = np.array(mesh.elements)
#            x-coords     y-coords connectivity
plt.triplot(points[:,0], points[:,1], tri)
plt.show()
