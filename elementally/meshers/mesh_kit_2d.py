# This is the module file for meshing in 2D.
# It will contain a number of functions creating
# specific meshs, as well as other function
# creating subfeatures of boundaries; e.g.
# functions creating line segments or circle arcs.
#
# Created 29.6 2014.

import meshpy.triangle as triangle
import numpy as np
import numpy.linalg as la

#---------------------------------------------|
#   Functions creating specific meshs.        |
#---------------------------------------------|

def quarter_annulus_2d(volume_tolerance, angle_start, angle_end, center, radius_inner, radius_outer):
    """
    Function creating a mesh for a quarter of an annulus
    lying on the first quadrant. For now, outer radius is set
    to 2.0 and inner radius 1.0.
    INPUT:
        volume_tolerance: Scalar specifying upper bound for area
            of each triangle in the mesh.
    OUTPUT:

    Typical usage:
    """

    # preprocessing
    num_points_outer = int( np.ceil(radius_outer * np.abs(angle_end - angle_start) /
                               np.sqrt(2 * volume_tolerance)) )
    num_points_inner = int( np.ceil(radius_inner * np.abs(angle_end - angle_start) /
                               np.sqrt(2 * volume_tolerance)) )

    facet_markers = []
    # Create points for outer radius:
    points = circle_segment(angle_start, angle_end, center, radius_outer, num_points_outer)
    # Boundary markers:
    #   1: Neumann boundary.
    #   2: Dirichlet boundary.
    facet_markers.extend( [1] * (num_points_outer - 1) )
    # Extending points list to include inner radius:
    points.extend( circle_segment(angle_end, angle_start, center, radius_inner, num_points_inner) )
    facet_markers.append(2)
    facet_markers.extend( [1] * (num_points_inner - 1) )
    facet_markers.append(2)
    num_points = len(points)
    # Create list of point connectivity:
    facets = connect_points(0,num_points-1)
    # Connect end point to starting point:
    facets.append( (num_points-1, 0) )
    
    #Use meshpy.triangle to create mesh:
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets, facet_markers)
    return triangle.build(info, max_volume = volume_tolerance)


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
    return [(i,i+1) for i in range(start_index, end_index)]

def facet_interior_point(facet, element_array):
    """
    For a vector facet on the form of (p0, p1), return the remaining third index (the interior point)
    from the element array of the mesh.
    INPUT:
        facet - tuple of two indeces
        element_array - array (n, 3) size, where n is the number of elements in the mesh
    OUTPUT:
        integer value of remaining index
    """
    for element in element_array:
        if np.all(np.intersect1d(facet, element) == np.sort(facet)):
            return np.setdiff1d(element, facet)[0]
    else:
        print "Did not find any elements containing both indices {0}".format(facet)

def check_facet_direction(facet, element_array, points):
    """
    For a vector facet (i,j), check if going from point i to point j,
    is in the counterclock direction.
    INPUT:
        facet - tuple of two indices.
        element_array - array (n, 3) size, where n is the number of elements in the mesh.
        points - array (m, 2) size, where m is the number of points in the mesh.
    OUTPUT:
        Returns True if going from i to j is in the counterclockwise direction.
    """
    point1 = points[facet[0],:]
    point2 = points[facet[1],:]
    point3 = points[facet_interior_point(facet,element_array),:]
    T = np.vstack ( (point2-point1, point3-point1) )
    return ( la.det(T) > 0.)

def facet_orientation(facets, element_array, points):
    """
    Checks direction on all facets, and returns an updated list
    of facets.
    INPUT:
        facets - list of tuples with facets.
        elements_array - n*3 array, where n is number of elements.
        points - m*2 array, where m is number of points.
    OUTPUT:
        updated_facets - list of facets with updated orientation.
    """
    updated_facets = []
    for facet in facets:
        if check_facet_direction(facet, element_array, points):
            updated_facets.append(facet)
        else:
            updated_facets.append( [facet[1], facet[0] ] )
    return updated_facets
