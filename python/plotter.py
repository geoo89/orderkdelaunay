import itertools
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import numpy as np
import pylab

class Plotter:
    '''Abstract class for plotting order-k Delaunay mosaics.

    The actual plotting functionality is implemented in Plotter2D and
    Plotter3D. All classes use the constructor from Plotter, which
    defines the following public attributes for plotting settings:

    Public attributes:
        draw_input_points: True/False - Whether to draw the input point set.
        draw_labels: True/False - Whether to draw labels for the vertices and
                     input points.
        colors_cells: list of 3 matplotlib.colors - for first, second and
                      third generation cells.
    '''

    def __init__(self, points, orderk_delaunay):
        '''Initialize various drawing settings.

        These are public attributes and can be modified by the user.
        '''
        self._points = points
        self._orderk_delaunay = orderk_delaunay

        # Whether to draw a label for each vertex
        self.draw_labels = True

        # Whether to draw a label for each vertex
        self.draw_input_points = True

        # colors of the cells of each generation
        color_firstgen_cell = colors.rgb2hex((1.0, 0.0, 0.0))
        color_secondgen_cell = colors.rgb2hex((0.0, 0.0, 1.0))
        color_thirdgen_cell = colors.rgb2hex((1.0, 1.0, 0.0))
        self.colors_cells = [
              color_firstgen_cell, color_secondgen_cell, color_thirdgen_cell]

        # separator of the points indices that the vertex label consists of
        self._label_sep = ''
        if len(points) > 10:
            self._label_sep = ','

    def draw(self, order):
        '''Stub. Implemented by subclasses.'''
        pass


class Plotter2D(Plotter):

    def __init__(self, points, orderk_delaunay):
        super().__init__(points, orderk_delaunay)

    def draw(self, order):
        vertices = self._orderk_delaunay.diagrams_vertices[order-1]
        cells = self._orderk_delaunay.diagrams_cells[order-1]
        generations = self._orderk_delaunay.diagrams_generations[order-1]

        ax = plt.gca()
        ax.cla()
        # Draw cells stemming from barycentric subdivisions of older cells.
        # Either new_nextgen_cells or triangulated_cells is empty, depending on
        # whether 'triangulate' is True or False.
        for tri, gen in zip(cells, generations):
            vxs = []
            for vx in tri:
                xs = [self._points[c][0] for c in vx]
                ys = [self._points[c][1] for c in vx]
                xcenter = sum(xs)/len(xs)
                ycenter = sum(ys)/len(ys)
                vxs.append([xcenter, ycenter])
            p = plt.Polygon(vxs, closed=True, fill=True,
                            color=self.colors_cells[gen-1], alpha = 0.4)
            ax.add_patch(p)
            for pair in itertools.combinations(vxs, 2):
                pair = np.array(pair)
                line = plt.Line2D(pair[:, 0], pair[:, 1], 
                      color=self.colors_cells[gen-1], linewidth = 1.0)
                ax.add_artist(line)

        # order-k points
        centroids = np.array([np.sum(
              [self._points[i] for i in vertex], axis=0)/order
                    for vertex in vertices])
        ax.plot(centroids[:,0], centroids[:,1], 'o', color="black")
        if self.draw_labels:
            for c, v in zip(centroids, vertices):
                ax.text(c[0], c[1], self._label_sep.join([str(pt) for pt in v]))
        # first-order points
        if self.draw_input_points:
            ax.plot(self._points[:,0], self._points[:,1], 'o', color='red')
            if self.draw_labels:
                for i in range(len(self._points)):
                    ax.text(self._points[i][0], self._points[i][1], str(i))

        plt.show()


class Plotter3D(Plotter):

    def __init__(self, points, orderk_delaunay):
        super().__init__(points, orderk_delaunay)

        # Factor by how much to shrink simplices (helps visualization).
        self.shrinking_factor = 0.25

    def draw(self, order):
        vertices = self._orderk_delaunay.diagrams_vertices[order-1]
        cells = self._orderk_delaunay.diagrams_cells[order-1]
        generations = self._orderk_delaunay.diagrams_generations[order-1]

        ax = a3.Axes3D(pylab.figure())
        for cell, gen in zip(cells, generations):
            gvertices = []
            # compute the geometric vertices
            for gvertex in cell:
                gvertices.append(np.mean(np.array(
                      [self._points[i][0:3] for i in gvertex]), axis=0))
            gvertices = np.array(gvertices)
            col = self.colors_cells[gen-1]
            faces = []
            cell_center = np.mean(gvertices, axis=0)

            # tetrahedron
            if len(cell) == 4:
                for facet in itertools.combinations(gvertices, 3):
                    # geometric facet
                    gfacet = [gv[0:3]*(1-self.shrinking_factor) + \
                          cell_center*(self.shrinking_factor) for gv in facet]
                    faces.append(gfacet)
            # octahedron
            elif len(cell) == 6:
                for facet_ids in itertools.combinations(range(len(gvertices)), 3):
                    facet = [cell[i] for i in facet_ids]
                    gfacet = [gvertices[i] for i in facet_ids]
                    is_facet = True
                    # three vertices only span a facet if they pairwise differ
                    # in only one element.
                    for edge in itertools.combinations(facet, 2):
                        if len(set(edge[0]) & set(edge[1])) != order-1:
                            is_facet = False
                    if is_facet:
                        # geometric facet
                        gfacet = [gv[0:3]*(1-self.shrinking_factor) + \
                              cell_center*(self.shrinking_factor)
                                    for gv in gfacet]
                        faces.append(gfacet)
            tri = a3.art3d.Poly3DCollection(faces, alpha=0.2)
            tri.set_color(col)
            tri.set_edgecolor('k')
            ax.add_collection3d(tri)
        # order-k points
        centroids = np.array([np.sum(
              [self._points[i] for i in vertex], axis=0)/order
                    for vertex in vertices])
        ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2], 'o',
                   color="black")
        if self.draw_labels:
            for c, v in zip(centroids, vertices):
                ax.text(c[0], c[1], c[2], self._label_sep.join([str(pt) for pt in v]))
        # first-order points
        if self.draw_input_points:
            ax.scatter(self._points[:,0], self._points[:,1], self._points[:,2],
                  marker='o',color='red')
            if self.draw_labels:
                for i in range(len(self._points)):
                    ax.text(self._points[i][0], self._points[i][1],
                            self._points[i][2], str(i))

        plt.show()
