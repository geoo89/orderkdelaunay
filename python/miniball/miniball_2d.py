import math
import numpy as np

# numerical error margin
EPS = 1e-12


def circumsphere_2d(points):
    '''Circumsphere of up to 3 points.

    Args:
        A list of the input points.

    Returns:
        Pair of circumcenter and circumradius.
    '''

    if len(points) == 0:
        return None, 0
    elif len(points) == 1:
        return points[0], 0
    elif len(points) == 2:
        return (np.array(points[0]) + np.array(points[1])) / 2, \
                0.5*np.linalg.norm(np.array(points[0]) - np.array(points[1]))
    elif len(points) == 3:
        # Mathematical algorithm from Wikipedia: Circumscribed circle
        ax = points[0][0]; ay = points[0][1]
        bx = points[1][0]; by = points[1][1]
        cx = points[2][0]; cy = points[2][1]
        d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
        if d == 0.0:
            return None
        x = ((ax * ax + ay * ay) * (by - cy) + \
             (bx * bx + by * by) * (cy - ay) + \
             (cx * cx + cy * cy) * (ay - by)) / d
        y = ((ax * ax + ay * ay) * (cx - bx) + \
             (bx * bx + by * by) * (ax - cx) + \
             (cx * cx + cy * cy) * (bx - ax)) / d
        return (x, y), math.hypot(x - ax, y - ay)
    else:
        raise ValueError("Can't have 2D-circumsphere of more than 3 points.")


def miniexball_2d(pin, pon, pout):
    """Compute smallest enclosing ball of pin that has pon
    on its surface and pout outside or on the surface.

    Args:
        pin: List of points to be inside the ball
        pon: List of points to be on the ball
        pout: List of points to be outside the ball

    Returns:
        tuple: the center (x, y) of the enclosing ball
        int: the radius of the enclosing ball
    """

    if len(pon) == 3:
        cc, cr = circumsphere_2d(pon)
        valid = True
        for p in pin:
            if np.linalg.norm(np.array(p) - np.array(cc)) > cr + EPS:
                valid = False
        for p in pout:
            if np.linalg.norm(np.array(p) - np.array(cc)) < cr - EPS:
                valid = False
        if not valid:
            raise ValueError("No sphere including pin and excluding pout exists.")
    elif pout != []:
        p = pout.pop()
        # compute smallest enclosing disk of pin-p, pon
        cc, cr = miniexball_2d(list(pin), list(pon), list(pout))
        # if p is inside the disk
        if cc is None or np.linalg.norm(np.array(p) - np.array(cc)) < cr:
            cc, cr = miniexball_2d(list(pin), list(pon + [p]), list(pout))
    elif pin != []:
        p = pin.pop()
        # compute smallest enclosing disk of pin-p, pon
        cc, cr = miniexball_2d(list(pin), list(pon), list(pout))
        # if p is outside the disk
        if cc is None or np.linalg.norm(np.array(p) - np.array(cc)) > cr:
            cc, cr = miniexball_2d(list(pin), list(pon + [p]), list(pout))
    else:
        cc, cr = circumsphere_2d(pon)
    return cc, cr
