import itertools
import numpy as np
import phat
import scipy.spatial
from miniball import miniball
from orderk_delaunay import OrderKDelaunay

'''
Algorithm to compute persistence of the k-fold cover for a given point set.

Outline:
- Compute order-k Delaunay triangulation (top-dimensional simplices)
- Make dict of all simplices (any dimension), containing:
      co-faces and dimension
- Compute filtration value for all simplices using miniball variant
      (This needs the coface information)
- Sort somplices by dimension and filtration value
- Assign indices to each simplex and make boundary matrix
- invoke phat to compute persistence
'''


class CellInfo:
    '''
    Struct to store information for a combinatorial cell.

    Attributes:
        coface_vxs: Vertices that can be added to the cell to form a coface.
        radius: radius/filtation value of the cell.
        dimension: dimension of the cell.
        bdmx_index: Index of the cell in the boundary matrix.
    '''

    def __init__(self, coface_vxs, radius, dimension, bdmx_index):
        self.coface_vxs = coface_vxs
        self.radius = radius
        self.dimension = dimension
        self.bdmx_index = bdmx_index


def kcover_persistence(points, order):
    '''
    Compute persistence of the k-fold cover of balls for a set of points
    in 2D or 3D.

    Args:
        points: list of points
        order: order k of the k-fold cover with respect to which 
               to compute persistence

    Returns:
        ppairs:
            List of persistence pairs.
        filtration:
            Dictionary mapping each simplex of the triangulated order-k
            Delaunay triangulation to a CellInfo struct.
        filtration_sorted:
            Same as filtration, but sorted by dimension with filtration
            value as tiebreaker.
    '''


    # Dimension of the ambient space.
    dimension = len(points[0])

    # Compute order-k Delaunay triangulation
    okdel = OrderKDelaunay(points, order)
    ktuples_allk = okdel.diagrams_vertices
    simplices_allk = okdel.diagrams_simplices
    # Vertex set (represented as k-tuples of input points)
    ktuples = ktuples_allk[order-1]
    # Top-dimensional simplices.
    # Simplices are tuples of indices into the vertex set.
    simplices = simplices_allk[order-1]


    # From top-dimensional cells in decreasing dimension down to vertices,
    # make a dictionary mapping each cell (as tuple of its vertices) to
    # a CellInfo instance containing its co-face vertices, filtration value,
    # dimension and index in the boundary matrix.
    filtration_topdimcells = dict()
    simplices = [tuple(simplex) for simplex in simplices]
    for simplex in simplices:
        filtration_topdimcells[simplex] = CellInfo([], -1, dimension, -1)

    if dimension == 2:
        filtration_faces = filtration_topdimcells

        # We don't have tetrahedra in 2D, so this stays empty.
        filtration_tetras = dict()

    elif dimension == 3:
        filtration_tetras = filtration_topdimcells

        filtration_faces = dict()
        for tetr in filtration_tetras:
            # Selection of three vertices.
            # The first three entries are the vertices making the face. The last
            # one is the vertex together with which the face makes a co-face
            for sel in [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [1, 2, 3, 0]]:
                # get the face consisting of three vertices
                face = (tetr[sel[0]], tetr[sel[1]], tetr[sel[2]])
                # add the fourth (missing) vertex to its dictionary entry
                # together with this vertex we get a co-face
                # (which is our tetr that we started with)
                if face in filtration_faces.keys():
                    filtration_faces[face].coface_vxs.append(tetr[sel[3]])
                else:
                    filtration_faces[face] = CellInfo([tetr[sel[3]]], -1, 2, -1)

    filtration_edges = dict()
    for face in filtration_faces:
        # Selection of two vertices.
        for sel in [[0, 1, 2], [0, 2, 1], [1, 2, 0]]:
            # get the edge consisting of two vertices
            edge = (face[sel[0]], face[sel[1]])
            # add the third (missing) vertex to its dictionary entry
            # together with this vertex we get a co-face
            # (which is our face that we started with)
            if edge in filtration_edges.keys():
                filtration_edges[edge].coface_vxs.append(face[sel[2]])
            else:
                filtration_edges[edge] = CellInfo([face[sel[2]]], -1, 1, -1)

    filtration_vxs = dict()
    for edge in filtration_edges:
        for sel in [[0, 1], [1, 0]]:
            vertex = (edge[sel[0]], )
            if vertex in filtration_vxs:
                filtration_vxs[vertex].coface_vxs.append(edge[sel[1]])
            else:
                filtration_vxs[vertex] = CellInfo([edge[sel[1]]], -1, 0, -1)

    # Merge cells of all dimensions into a single new dict.
    filtration = dict()
    filtration.update(filtration_tetras)
    filtration.update(filtration_faces)
    filtration.update(filtration_edges)
    filtration.update(filtration_vxs)


    # Compute the radius values for each simplex
    for simplex in filtration:
        # flatten the simplex (i.e. the list of ktuples of points):
        # We just take all the ktuples and put them in a single list.
        flattened = itertools.chain.from_iterable([ktuples[i] for i in simplex])
        # set of points in or on the sphere
        pinon = set(flattened)
        # set of points in (possibly on) the sphere are those
        # that appear in every vertex of our simplex
        pin = pinon.intersection(*[ktuples[i] for i in simplex])
        # Those that don't appear in every vertex have to be ON the sphere,
        # not inside.
        pon = pinon - pin
        # flatten the list of vertices (i.e. point ktuples) that are
        # part of co-faces of our simplex
        flattened = itertools.chain.from_iterable(
              [ktuples[i] for i in filtration[simplex].coface_vxs])
        # These are the points that must be outside (or on) the sphere
        pout = set(flattened) - pinon

        cc, cr = miniball.miniexonball(
              [points[i] for i in pin], 
              [points[i] for i in pon], 
              [points[i] for i in pout])
        # Assign the filtration value.
        filtration[simplex].radius = cr


    # Sort by dimension, then by filtration value as tiebreaker.
    # This ensures a cell will always appear after its faces,
    # even in case of rounding errors.
    filtration_sorted = sorted(filtration.items(),
          key=lambda x: (x[1].dimension, x[1].radius))

    # Assign each simplex its index in the list.
    # Note: Because the list of data appears both in the dictionary and
    # the sorted version as a pointers, changing its content in the sorted
    # version is sufficient to also change the content in the dict
    for i in range(len(filtration_sorted)):
        filtration_sorted[i][1].bdmx_index = i

    # Make boundary matrix of the complex.
    boundary_matrix = phat.boundary_matrix(
          representation = phat.representations.vector_vector)
    bdmx = []
    for simplex in filtration_sorted:
        dim = simplex[1].dimension
        if dim == 0:
            # If it's a vertex, it doesn't have a boundary
            bd = []
        else:
            # For each face, select all but 1 point.
            # Look up the index in the filtration dictionary and 
            # put the list of these indices as boundary.
            bd = [filtration[i].bdmx_index for i in 
                  itertools.combinations(simplex[0], dim)]
        bdmx.append((dim, sorted(bd)))
    boundary_matrix.columns = bdmx

    # Compute persistence using phat
    ppairs = boundary_matrix.compute_persistence_pairs()
    ppairs.sort()

    return ppairs, filtration, filtration_sorted
