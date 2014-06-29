# This is the module file for meshing in 2D.
# It will contain a number of functions creating
# specific meshs, as well as other function
# creating subfeatures of boundaries; e.g.
# functions creating line segments or circle arcs.
#
# Created 29.6 2014.

import meshpy.triangle as triangle
import numpy as np

#---------------------------------------------|
#   Functions creating specific meshs.        |
#---------------------------------------------|

def quarter_annulus_2D()
    """
    Function creating a mesh for a quarter of an annulus
    lying on the first quadrant.
    INPUT:

    OUTPUT:

    Typical usage:
    """
    pass




#---------------------------------------------|
#   Functions useful for creating mesh.       |
#---------------------------------------------|

def circle_segment(theta_start, theta_end, center, radius, num_points):
    """
    Function creating points lying on a circle arc.
    INPUT:
        theta_start: Start angle in radians.
        theta_end: End angle in radians.
        center: Point (x,y) specifying the center of the circle
                which the segment lies on.
        radius: Scalar specifying the radius of that circle.
        num_points: integer specifying the number of points
                    to be returned.
    OUTPUT:
        points: A list of uniformly distributed points
                along the circle arc. The points will be ordered
                according to angle. From theta_start to theta_end.

    Typical usage:
        circle_segment(0, np.pi/2., (0,0), 1.0, 10) returns
        10 points lying on the unit circle in the first quadrant,
        ordered in inreasing angle.
    """
    theta = np.linspace(theta_start, theta_end, num_points)
    # Initialise list of points:
    points = []
    # Iterating over each angle:
    for i in range(num_points):
        #              (        x-coordinate                 ,          y-coordinate                )
        points.append( (center[0] + radius * np.cos(theta[i]), center[1] + radius * np.sin(theta[i])) )
    # NOTE: This could be done in a single line using list comprehension (<3).
    return points

def connect_points(start_index, end_index):
    """
    Function returning list of point connectivity.
    INPUT:
        start_index and end_index specifying which points to be connected.
    OUTPUT:
        List of point connectivity.
    Typical usage:
        connect_points(0,3) would return the list
        [(0,1), (1,2), (2,3)].
        Meaning that point_0 is connected to point_1,
        point_1 connected to point_2 and so on.
    """
    return [(i,i+1) for i in range(start_index, end_index)
