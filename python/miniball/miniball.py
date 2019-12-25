from . import miniball_2d
from . import miniball_3d

"""
Functions to compute the smallest enclosing ball of a given set of
points, with possible additional constrains of points to be excluded
or to be on the surface of the ball.

This is an adaptation of the algorithm by Welzl from the 1991 paper
"Smallest enclosing disks (balls and ellipsoids)". Thanks to Grzegorz
Jabłoński for making me aware of its generalizability.

TODO: Factor out the dimension-dependent aspects of miniexball_2d
and miniexball_3d as the main logic is the same for 2D and 3D.
"""

def miniball(pin):
    """Compute smallest enclosing ball of pin in 2 or 3 dimensions.

    Args:
        pin: List of points to be inside/on the ball

    Returns:
        tuple: the center (x, y) of the enclosing ball
        int: the radius of the enclosing ball
    """
    if len(pin) != 0:
        dim = len(pin[0])
    else:
        return None, 0

    if dim == 2:
        return miniball_2d.miniexball_2d(pin, [], [])
    elif dim == 3:
        return miniball_3d.miniexball_3d(pin, [], [])
    else:
        raise ValueError("Miniball is not implemented for dimension %d" % dim)


def miniexball(pin, pout):
    """Compute smallest enclosing ball of pin which does not
       contain any points from pout in its interior.

    Args:
        pin: List of points to be inside/on the ball
        pout: List of points to be outside/on the ball

    Returns:
        tuple: the center (x, y) of the enclosing ball
        int: the radius of the enclosing ball

    Raises:
        ValueError: if no such ball exists
    """

    if len(pin) != 0:
        dim = len(pin[0])
    else:
        return None, 0

    if dim == 2:
        return miniball_2d.miniexball_2d(pin, [], pout)
    elif dim == 3:
        return miniball_3d.miniexball_3d(pin, [], pout)
    else:
        raise ValueError("Miniball is not implemented for dimension %d" % dim)


def miniexonball(pin, pon, pout):
    """Compute smallest enclosing ball of pin which does not
       contain any points from pout in its interior and has the
       points of pon on its boundary.

    Args:
        pin: List of points to be inside/on the ball
        pon: List of points to be on the ball
        pout: List of points to be outside/on the ball

    Returns:
        tuple: the center (x, y) of the enclosing ball
        int: the radius of the enclosing ball

    Raises:
        ValueError: if no such ball exists
    """

    if len(pin) != 0:
        dim = len(pin[0])
    elif len(pon) != 0:
        dim = len(pon[0])
    else:
        return None, 0
    # TODO: To deal with degeneracies, we should use the largest
    # affinely independent subset of pon rather than pon itself
    if dim == 2:
        return miniball_2d.miniexball_2d(pin, pon, pout)
    elif dim == 3:
        return miniball_3d.miniexball_3d(pin, pon, pout)
    else:
        raise ValueError("Miniball is not implemented for dimension %d" % dim)
