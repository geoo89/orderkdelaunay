import math
import numpy as np


# numerical error margin
EPS = 1e-12


def normal(points):
    '''Compute normal vector of two input vectors.

    Args:
        points: List of two lists of length 3
                representing two 3-dimensional points.
    
    Return:
        Normal vector of the two input points (vectors).
    '''

    p = points[0]
    q = points[1]

    n = np.array([p[1]*q[2] - p[2]*q[1],
                  p[2]*q[0] - p[0]*q[2],
                  p[0]*q[1] - p[1]*q[0]])
    return n


def circumsphere_3d_2(pointsf):
    '''Circumsphere of 2 points.

    Args:
        A list of the 2 input points.

    Returns:
        Pair of circumcenter and circumradius.
    '''

    points = np.array(pointsf)
    return list((points[0] + points[1]) / 2.0), \
           np.linalg.norm((points[0] - points[1]) / 2.0)


def circumsphere_3d_3(pointsf):
    '''Circumsphere of 3 points.

    Args:
        A list of the 3 input points.

    Returns:
        Pair of circumcenter and circumradius.
    '''

    points = np.array(pointsf)

    d1 = points[1] - points[0]
    d2 = points[2] - points[0]

    # compute normal vector
    n = normal([d1,d2])

    # Make system of 3 equations to solve for.
    # Third equation is: Circumcenter should be orthogonal to normal.
    genmat = np.array([d1, d2, n])
    normvec = np.array([np.linalg.norm(genmat[0])**2, 
                        np.linalg.norm(genmat[1])**2, 
                        0])
    # 0-based circumsphere is the solution to this equation system.
    cc = np.linalg.solve(2*genmat, normvec)

    return list(cc + points[0]), np.linalg.norm(cc)


def circumsphere_3d_4(pointsf):
    '''Circumsphere of 4 points.

    Args:
        A list of the 4 input points.

    Returns:
        Pair of circumcenter and circumradius.
    '''

    points = np.array(pointsf)

    # generator matrix D: Contains (as rows) the three direction vectors of
    # the tetrahedron spanned by the points with the first point as origin.
    genmat = np.array([points[i] - points[0] for i in range(1,4)])

    # vector n containing the square norms of the 3 vectors from genmat.
    normvec = np.array([np.linalg.norm(genmat[i])**2 for i in range(3)])

    # 0-based circumsphere is the solution x to 2Dx = n
    cc = np.linalg.solve(2*genmat, normvec)

    # Return translated circumcenter and 
    return list(cc + points[0]), np.linalg.norm(cc)


def circumsphere_3d(points):
    npoints = len(points)
    if npoints == 0:
        return None, 0
    elif npoints == 1:
        return points[0], 0
    elif npoints == 2:
        return circumsphere_3d_2(points)
    elif npoints == 3:
        return circumsphere_3d_3(points)
    elif npoints == 4:
        return circumsphere_3d_4(points)
    else:
        raise ValueError("Cannot have {} points on a sphere.".format(npoints))


def miniexball_3d(pin, pon, pout):
    """Compute smallest enclosing ball of pin that has pon
    on its surface.

    Args:
        pin: List of points to be inside the ball
        pon: List of points to be on the ball
        pout: List of points to be outside the ball

    Returns:
        tuple: the center (x, y, z) of the enclosing ball
        int: the radius of the enclosing ball
    """

    if len(pon) == 4:
        cc, cr = circumsphere_3d(pon)
        valid = True
        for p in pin:
            if np.linalg.norm(np.array(p) - np.array(cc)) > cr + EPS:
                valid = False
        for p in pout:
            if np.linalg.norm(np.array(p) - np.array(cc)) < cr - EPS:
                valid = False
        if not valid:
            raise ValueError("No sphere including {} and excluding {} exists."
                    .format(pin, pout))
    elif pout != []:
        p = pout.pop()
        # compute smallest enclosing disk of pin-p, pon
        cc, cr = miniexball_3d(list(pin), list(pon), list(pout))
        # if p is inside the disk
        if cc is None or np.linalg.norm(np.array(p) - np.array(cc)) < cr:
            cc, cr = miniexball_3d(list(pin), list(pon + [p]), list(pout))
    elif pin != []:
        p = pin.pop()
        # compute smallest enclosing disk of pin-p, pon
        cc, cr = miniexball_3d(list(pin), list(pon), list(pout))
        # if p is outside the disk
        if cc is None or np.linalg.norm(np.array(p) - np.array(cc)) > cr:
            cc, cr = miniexball_3d(list(pin), list(pon + [p]), list(pout))
    else:
        cc, cr = circumsphere_3d(pon)
    return cc, cr

