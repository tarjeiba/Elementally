# A short script showing how to use the mesh_kit_2D to create a mesh of quarter_annulus,
# and plots the resulting mesh using matplotlib.
#
# Created 29.6 2014

import numpy as np
import matplotlib.pyplot as plt
import mesh_kit_2d

mesh = mesh_kit_2d.quarter_annulus_2d(.1, np.pi/2, np.pi/4, (0, 0), 1.0, 2.0)
points = np.array(mesh.points)
tri = np.array(mesh.elements)
edges = np.array(mesh.facets)
print " Characteristics of mesh:"
print "  Number of nodes: " + str(points.shape[0])
print "  Number of elements: " + str(tri.shape[0])
print "  Number of boundary edges: " + str( edges.shape[0])
#            x-coords     y-coords connectivity
plt.triplot(points[:,0], points[:,1], tri)
plt.show()
print mesh.points