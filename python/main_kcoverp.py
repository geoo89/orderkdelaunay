import collections
import kcover_persistence as kcoverp
import matplotlib.pyplot as plt
import numpy as np

'''
Print persistence pairs and draw persistence diagram of the k-fold cover
for a given 2- or 3-dimensional point set

Colors of persistence pairs in persistence diagram:
Blue: 0-dimensional
Magenta: 1-dimensional
Red: 2-dimensional
'''

if __name__ == "__main__":
    # # 2D example point set
    # points = np.array([[156.006,705.854],
    #         [215.257,732.63], [283.108,707.272], [244.042,670.948],
    #         [366.035,687.396], [331.768,625.715], [337.936,559.92],
    #         [249.525,582.537], [187.638,556.13], [165.912,631.197]])

    # 3D example point set
    points = np.array([(0,0,0), (0,4,4), (4,4,0), (4,0,4), (-10,2,2)])

    # the order k of the k-fold cover
    order = 2

    # Rounding of the birth and death values of the persistence pairs.
    # This is useful to group the persistence pairs and assign multiplicities
    rounding_digits = 3
    # Draw the multiplicities of each point in the persistence diagram.
    draw_multiplicities = True
    # We only draw pairs of non-zero persistence. 
    # This is the cutoff above which we consider the persistence to be non-zero
    persistence_cutoff = 1e-3
    # How much to perturb the points. May be 0, but the code cannot deal with
    # degeneracies (i.e. points not in general position)
    # This is a uniform perturbation rather than Gaussian.
    perturbation_level = 1e-4

    dimension = len(points[0])
    # Perturb the points a little to avoid degeneracies
    if perturbation_level > 0:
        points = np.array(points) + np.random.random(np.array(points).shape) \
                 * perturbation_level - perturbation_level/2


    ppairs, _, filtration_sorted = kcoverp.kcover_persistence(points, order)

    pdict = collections.defaultdict(int)
    print("There are %d persistence pairs." % len(ppairs))
    print("The following ones have non-zero persistence.")
    for pair in ppairs:
        dim = filtration_sorted[pair[0]][1].dimension
        birth = round(filtration_sorted[pair[0]][1].radius, rounding_digits)
        death = round(filtration_sorted[pair[1]][1].radius, rounding_digits)
        if death - birth > persistence_cutoff:
            pdict[(birth, death, dim)] += 1

    for p,m in sorted(pdict.items()):
        print("Birth: %f, Death: %f, Dimension: %d, Multiplicity: %d" % \
              (p[0], p[1], p[2], m))

    # smallest and biggest birth/death time
    lower = min(pdict.keys(), key=lambda x: x[0])[0]
    upper = max(pdict.keys(), key=lambda x: x[1])[1]
    prange = upper - lower

    # draw the diagram
    ax = plt.gca()
    ax.cla()
    ax.set_aspect('equal', adjustable='box')
    ax.set_ylim((lower - 0.05*prange, upper + 0.05*prange))
    ax.set_xlim((lower - 0.05*prange, upper + 0.05*prange))
    for p,m in pdict.items():
        birth = p[0]
        death = p[1]
        dim = p[2]
        if death - birth > 1e-3:
            if dim == 0:
                color = 'blue'
            elif dim == 1:
                color = 'magenta'
            else:
                color = 'red'
            ax.plot([birth], [death], 'o', color=color)
            if draw_multiplicities:
                plt.text(birth, death, m,
                      color="black", fontsize=10, horizontalalignment='right')
    line = plt.Line2D([0,999999], [0,999999], color='black')
    ax.add_artist(line)
    plt.show()
