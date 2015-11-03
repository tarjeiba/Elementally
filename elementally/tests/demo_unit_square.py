#!/usr/bin/env python
#
# Demo of creating a mesh for the unit square
# and testing the all_edges feature, which
# will be important when considering mixed
# formulations of e.g. Poisson, Stokes and 
# Maxwell's equation.
##--------------------------------------##

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

import meshers

############################################

# Create mesh:
mesh = meshers.unit_square_2d(2,2, generate_facets=True, \
      generate_neighbors=True)

points = np.array(mesh.points)
faces = np.array(mesh.faces)
face_markers = np.array(mesh.face_markers)
print "Points: ", points
print "Faces: ", faces
print "Face markers: ", face_markers
print "Elements: ", np.array(mesh.elements)
print "Facets: ", np.array(mesh.facets)
print "Facet markers: ", np.array(mesh.facet_markers)

print "Neighbors: ", np.array(mesh.neighbors)
print "Normals: ", np.array(mesh.normals)

# Plot mesh:
fig = plt.figure(1)
ax = fig.gca()

ax.triplot(points[:,0], points[:,1],
          triangles=np.array(mesh.elements))

plt.show(1)
