import os
import sys
# Allow importing any modules relative to the main path.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
import numpy as np
from math import sqrt

from miniball.miniball_2d import circumsphere_2d
from miniball.miniball_2d import miniexball_2d
from miniball.miniball_3d import circumsphere_3d
from miniball.miniball_3d import miniexball_3d
from miniball.miniball import miniball, miniexball, miniexonball

class TestMiniball2d(unittest.TestCase):

    def test_cc2_2d(self):
        cc, r = circumsphere_2d([[0, 0], [2, 0]])
        np.testing.assert_allclose(cc, [1, 0])
        self.assertAlmostEqual(r, 1)

    def test_cc3_2d(self):
        cc, r = circumsphere_2d([[0, 0], [1, 0], [0, 1]])
        np.testing.assert_allclose(cc, [0.5, 0.5])
        self.assertAlmostEqual(r, sqrt(2)/2)

    def test_cc2_3d(self):
        cc, r = circumsphere_3d([[0, 0, 0], [2, 0, 0]])
        np.testing.assert_allclose(cc, [1, 0, 0])
        self.assertAlmostEqual(r, 1)

    def test_cc3_3d(self):
        cc, r = circumsphere_3d([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        np.testing.assert_allclose(cc, [0.5, 0.5, 0])
        self.assertAlmostEqual(r, sqrt(2)/2)

    def test_cc4_3d(self):
        cc, r = circumsphere_3d([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        np.testing.assert_allclose(cc, [0.5, 0.5, 0.5])
        self.assertAlmostEqual(r, sqrt(3)/2)

    def test_mini1_2d(self):
        cc, r = miniball([[0, 0], [1, 0], [0.5, 0.5]])
        np.testing.assert_allclose(cc, [0.5, 0.0])
        self.assertAlmostEqual(r, 0.5)

    def test_mini2_2d(self):
        cc, r = miniball([[0, 0], [1, 0], [0.5, 1]])
        np.testing.assert_allclose(cc, [0.5, 0.375])
        self.assertAlmostEqual(r, 0.625)

    def test_miniex_2d(self):
        cc, r = miniexball([[0, 0], [1, 0]], [[0.5, 0.25]])
        np.testing.assert_allclose(cc, [0.5, -0.375])
        self.assertAlmostEqual(r, 0.625)

    def test_miniexon_2d(self):
        cc, r = miniexonball([[0, 0], [1, 0]], [[0.5, -1]], [[0.5, 0.3]])
        np.testing.assert_allclose(cc, [0.5, -0.375])
        self.assertAlmostEqual(r, 0.625)

    def test_mini_3d(self):
        cc, r = miniball([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        np.testing.assert_allclose(cc, [1/3, 1/3, 1/3])
        self.assertAlmostEqual(r, sqrt(2/3))

    def test_mini2_3d(self):
        cc, r = miniball([[-100, -100, -100], [55, 1, 0], [0, 33, 33], [100, 100, 100]])
        np.testing.assert_allclose(cc, [0, 0, 0])
        self.assertAlmostEqual(r, 100*sqrt(3))

    def test_miniex_3d(self):
        cc, r = miniexball([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]])
        np.testing.assert_allclose(cc, [0.5, 0.5, 0.5])
        self.assertAlmostEqual(r, sqrt(3)/2)

    def test_miniexon1_3d(self):
        cc, r = miniexonball([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]], [])
        np.testing.assert_allclose(cc, [0.5, 0.5, 0.5])
        self.assertAlmostEqual(r, sqrt(3)/2)

    def test_miniexon2_3d(self):
        cc, r = miniexonball([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]], [], [[0.5, 0.5, 0.5]])
        np.testing.assert_allclose(cc, [-0.25, -0.25, -0.25])
        self.assertAlmostEqual(r, sqrt(27)/4)

    def test_miniexon_planar_3d(self):
        cc, r = miniexonball([[0, 0, 0], [1, 0, 0]], [[0.5, -1, 0]], [[0.5, 0.3, 0]])
        np.testing.assert_allclose(cc, [0.5, -0.375, 0])
        self.assertAlmostEqual(r, 0.625)

    def test_miniexon_error_3d(self):
        with self.assertRaises(ValueError):
            miniexonball([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 0, 0]], [[0.5, 0.5, 0.5]])


if __name__ == '__main__':
    unittest.main()
