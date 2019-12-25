import numpy as np
from orderk_delaunay import OrderKDelaunay
from plotter import Plotter, Plotter2D, Plotter3D


if __name__ == "__main__":
    # # 2D example point set
    # points = np.array([[156.006,705.854],
    #         [215.257,732.63], [283.108,707.272], [244.042,670.948],
    #         [366.035,687.396], [331.768,625.715], [337.936,559.92],
    #         [249.525,582.537], [187.638,556.13], [165.912,631.197]])

    # 3D example point set
    points = np.array([(0,0,0), (0,4,4), (4,4,0), (4,0,4), (-10,2,2)])

    # the order k up to which to compute the order-k Delaunay diagram
    order = 4
    # Whether to print the cells of all the complexes
    print_output = True
    # Whether to draw all the order-k Delaunay mosaics
    draw_output = True

    # Compute the order-k Delaunay mosaics
    orderk_delaunay = OrderKDelaunay(points, order)

    # Initialize appropriate plotter for drawing the mosaics.
    dimension = len(points[0])
    if dimension == 2:
        plotter = Plotter2D(points, orderk_delaunay)
    elif dimension == 3:
        plotter = Plotter3D(points, orderk_delaunay)
    else:
        # Stub that doesn't draw anything.
        plotter = Plotter(points, orderk_delaunay)

    for k in range(1, order+1):
        if draw_output:
            plotter.draw(k)

        if print_output:
            cells = orderk_delaunay.diagrams_cells[k-1]
            # Output all the cells.
            print("Order {}. Number of cells: {}".format(
                  len(cells[0][0]), len(cells)))
            for cell in sorted(cells):
                print(cell)
