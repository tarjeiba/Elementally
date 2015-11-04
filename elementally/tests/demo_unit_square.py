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
mesh = meshers.unit_square_2d(2,2, generate_faces=True, \
      generate_neighbors=True)


points = np.array(mesh.points)
print "Points: ", points
print "Elements: ", np.array(mesh.elements)
print "Facets: ", np.array(mesh.facets)
print "Facet markers: ", np.array(mesh.facet_markers)

if mesh.faces:
    print "Faces: ", np.array(mesh.faces)
    print "Face markers: ", np.array(mesh.face_markers)
    #print "Neighbors: ", np.array(mesh.neighbors)
    #print "Normals: ", np.array(mesh.normals)

    # Let's try the new edges_opposite_vertex:
    mesh.set_edges_of_elements()
    print "Element edges: ", np.array(mesh.element_edges)

# Plot mesh:
fig = plt.figure(1)
ax = fig.gca()

ax.triplot(points[:,0], points[:,1],
          triangles=np.array(mesh.elements))

plt.show(1)
