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

from collections import Counter     # To be used in edge_opposite_vertex

class ElementallyMeshInfo(triangle.MeshInfo):

    # Attributes that can be set in this class.
    element_edges = []
    faces = []
    
    def facet_interior_point(self, facet):
        """
        For a vector facet on the form of (p0, p1), return the remaining third
        index (the interior point) from the element array of the mesh.
        INPUT:
            facet - tuple of two indeces
            element_array - array (n, 3) size, where n is the number of
            elements in the mesh
        OUTPUT:
            integer value of remaining index
        """
        for element in self.elements:
            # np.intersect1d: returns unique elements that are in both facet
            # and element.
            #
            # np.all both points in facet are element.
            #
            # np.setdiff1d give out index in element, not in facet.
            if np.all(np.intersect1d(facet, element) == np.sort(facet)):
                return np.setdiff1d(element, facet)[0]
        else:
            print "Did not find any elements containing both indices {0}"\
              .format(facet)

    def check_facet_direction(self, facet):
        """
        For a vector facet (i,j), check if going from point i to point j,
        is in the counterclock direction.
        INPUT:
            facet - tuple of two indices.
            element_array - array (n, 3) size, where n is the number of elements
            in the mesh.
            points - array (m, 2) size, m is the number of points in the mesh.
        OUTPUT:
            Returns True if going from i to j is in the counterclockwise
            direction.
        """
        point1 = np.array(self.points[facet[0]])
        point2 = np.array(self.points[facet[1]])
        point3 = np.array(self.points[(self.facet_interior_point(facet))])
        T = np.vstack ( (point2-point1, point3-point1) )
        return ( la.det(T) > 0.)

    def order_facets(self):
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

        for facet in self.facets:
            if self.check_facet_direction(facet):
                updated_facets.append(facet)
            else:
                updated_facets.append( [facet[1], facet[0] ] )

        # Something went wrong with using set_facets, so using __setstate__ instead.
        self.__setstate__((0,0,[["facets",updated_facets]]))
        #self.set_facets(updated_facets, self.facet_markers)


    def set_neighbors_from_voronoi(self, vorout):
        """
        Function for using a Voronoi diagram data structure to get
        neighbors and face normals of the mesh.
        INPUT:
          mesh: Output mesh from a call to build(ElementallyMeshInfo)
          vorout: Output voronoi diagram structure after call to Triangle.triangulate.
        OUTPUT:
          void: mesh is manipulated
        """
        if not (self.faces):
            print "mesh.faces is not initialized."
            return None

        neighbors = []
        normals = []
        for i, face in enumerate(self.faces):
            # Append neighboring triangles of this edge:
            neighbors.append((vorout.faces[i][0], vorout.faces[i][1]))      
            # Now to the normal of the edge to be pointing out of the triangle in
            # neighbors[i][0]:
            dirx = self.points[face[1]][0] - self.points[face[0]][0]
            diry = self.points[face[1]][1] - self.points[face[0]][1]
            mag = np.sqrt(dirx**2 + diry**2)

            normals.append( (diry/mag, -dirx/mag) ) 

        # Set states:
        self.__setstate__((0,0,[["neighbors", neighbors],\
                                ["normals", normals]]))
        #self.neighbors = neighbors
        #self.normals = normals

    def edge_opposite_vertex(self, vertex, element):
        """
        Function finding edge opposite vertex in element.
        Returns index of that edge.
        """
        if vertex not in self.elements[element]:
            print "Vertex is not in element."
            return -1

        # Get edge definition (end points):
        edge_def = [ind for ind in self.elements[element] if ind != vertex]
        # We need to run through all edges...
        for i, edge in enumerate(self.faces):
            # Check if the edges are the same:
            if (Counter(edge)==Counter(edge_def)):
                return i

        # Couldn't find edge:
        print "Couldn't find edge."
        return -1
    def set_edges_of_elements(self):
        """
        Function for generating an array that specifies
        what edges  constitute the boundary of an
        element. The result is ordered.
        OUTPUT:
          void, but element_edges is set. For instance
          element_edges[i] is a triple of integers where
          element_edges[i][0] is the global edge number
          opposite the vertex self.elements[i][0], and so on.
        """
        # Make sure that faces are set
        if not self.faces:
            print "Mesh faces is not initialized,"
            print "ElementallyMeshInfo::set_edges_of_elements() exits"
            print "without doing anything."
            return None

        # Iterate over all elements:
        element_edges = []
        for i, element in enumerate(self.elements):
            # For each element we run through each vertex,
            # and find opposite edge:
            temp = []
            for indv in xrange(3):
                vertex = element[indv]
                temp.append(self.edge_opposite_vertex(vertex, i))

            # Add this triple to element_edges:
            element_edges.append(temp)

        self.element_edges = element_edges

######################################################
##
##            BUILD METHOD
##
######################################################
def build(mesh_info, verbose=False, refinement_func=None, attributes=False,
        volume_constraints=False, max_volume=None, allow_boundary_steiner=True,
        allow_volume_steiner=True, quality_meshing=True,
        generate_edges=None, generate_faces=False, min_angle=None,
        mesh_order=None, generate_neighbors=False):
    """Triangulate the domain given in `mesh_info'.

    Taken with minute changes from triangle.build() written by Andreas Kloeckner.
    """
    opts = "pzj"
    if quality_meshing:
        if min_angle is not None:
            opts += "q%f" % min_angle
        else:
            opts += "q"

    if mesh_order is not None:
        opts += "o%d" % mesh_order

    if verbose:
        opts += "VV"
    else:
        opts += "Q"

    if attributes:
        opts += "A"

    if volume_constraints:
        opts += "a"
        if max_volume:
            raise ValueError, "cannot specify both volume_constraints and max_area"
    elif max_volume:
        opts += "a%.20f" % max_volume

    if refinement_func is not None:
        opts += "u"

    if generate_edges is not None:
        from warnings import warn
        warn("generate_edges is deprecated--use generate_faces instead")
        generate_faces = generate_edges

    if generate_faces:
        opts += "e"
        # THIS IS NEW: To get neighbors and normals
        if generate_neighbors:
            opts += "v" 
      

    if not allow_volume_steiner:
        opts += "YY"
        if allow_boundary_steiner:
            raise ValueError("cannot allow boundary Steiner points when volume "
                    "Steiner points are forbidden")
    else:
        if not allow_boundary_steiner:
            opts += "Y"

    # restore "C" locale--otherwise triangle might mis-parse stuff like "a0.01"
    try:
        import locale
    except ImportError:
        have_locale = False
    else:
        have_locale = True
        prev_num_locale = locale.getlocale(locale.LC_NUMERIC)
        locale.setlocale(locale.LC_NUMERIC, "C")

    try:
        mesh = ElementallyMeshInfo()
        vorout = triangle.MeshInfo()  # To be used to get neighbors
        # Interface for triangulate in backend Triangle is: 
        #   triangulate(options, input_info, output_mesh, voronoi_diagram)
        triangle.internals.triangulate(opts, mesh_info, mesh,
                                       vorout, refinement_func)

        # If generate_neighbors is true, we alse want neighbors and normals:
        if generate_neighbors and generate_faces:
            mesh.set_neighbors_from_voronoi(vorout)

    finally:
        # restore previous locale if we've changed it
        if have_locale:
            locale.setlocale(locale.LC_NUMERIC, prev_num_locale)

    return mesh

###############################################################################
##
##               GEOMETRIES AND BOUNDARIES
##
###############################################################################

def quarter_annulus_2d(volume_tolerance, angle_start, angle_end,
                       center, radius_inner, radius_outer):
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
    facet_markers.append(3)
    num_points = len(points)
    # Create list of point connectivity:
    facets = connect_points(0,num_points-1)
    # Connect end point to starting point:
    facets.append( (num_points-1, 0) )

    #Use meshpy.triangle to create mesh:
    info = ElementallyMeshInfo()
    info.set_points(points)
    info.set_facets(facets, facet_markers)
    mesh = build(info, max_volume = volume_tolerance)
    mesh.order_facets()
    return mesh

def unit_square_2d(nx, ny, generate_faces=False,
                  generate_neighbors=False):
  """
  Function creating a 2D mesh of the unit square [0,1]^2.
  INPUT:
    nx, ny: ints, specifying number of subdivisions in the
      x- and y direction, respectively.
    generate_facets: boolean specying whether all edges, not only boundary edges, should be
      generated in the output mesh.
  OUTPUT:
    mesh: Output mesh  
  """
  # Initialize point list:
  points = []
  
  # Add corner points:
  points.append( (0., 0.) )
  points.append( (1., 0.) )
  points.append( (1., 1.) )
  points.append( (0.,1.) )

  # Initialize segment list:
  facets = []
  facets = connect_points(0, 3)
  facets.append( (3, 0) )

  # Facet markers (for boundary):
  facet_markers = []
  facet_markers.extend( [1] * 4 )

  # Use meshpy.triangle to create mesh:
  info = ElementallyMeshInfo()
  info.set_points(points)
  info.set_facets(facets, facet_markers)
  mesh = build(info, max_volume = 0.5/(float(nx)*float(ny)),
        generate_faces=generate_faces,
        generate_neighbors=generate_neighbors)
  
  #Order boundary facets:
  mesh.order_facets()
  return mesh

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
