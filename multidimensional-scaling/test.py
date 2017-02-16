# -*- coding: utf-8 -*-
import unittest
import numpy as np
import scipy.io as sio
import sys
from numpy.linalg import eig

import utilities as util

if sys.version_info < (3, 3):
    from mock import MagicMock
else:
    from unittest.mock import MagicMock

class TestUtilities(unittest.TestCase):

    def setUp(self):
        # Loads data from any data file
        self.numberOfAnchors = 8
        self.anchors = np.random.randint(0,10,(self.numberOfAnchors,3))

    def test_distance_matrix(self):
        # Fixture
        A = np.array([1,2,3])
        B = np.array([4,2,6])
        C = np.array([2,3,9])
        P = np.array([A,B,C])
        dab = np.linalg.norm(A - B) ** 2
        dac = np.linalg.norm(A - C) ** 2
        dbc = np.linalg.norm(B - C) ** 2
        Dsq = np.array([[0,dab,dac],[dab,0,dbc],[dac,dbc,0]])
        
        # Test
        Dsq_comp = util.get_distance_matrix_squared(P)
        
        # Assert
        status = False
        try:
            np.testing.assert_array_almost_equal(Dsq, Dsq_comp, decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)

    def test_multi_dimensional_scaling(self):
        # Fixture
        Dsq = util.get_distance_matrix_squared(self.anchors)
        
        # Test
        Phat = util.multi_dimensional_scaling(Dsq)
        Dsq_comp = util.get_distance_matrix_squared(Phat)

        # Assert
        status = False
        try:
            np.testing.assert_array_almost_equal(Dsq, Dsq_comp, decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)
    
    def test_get_skew_symmetric_operator(self):
        # Fixture
        u = np.random.rand(3)
        v = np.random.rand(3)
        uv = np.cross(u,v)

        # Test
        S = util.get_skew_symmetric_operator(u)
        uv_comp = S.dot(v)
        
        # Assert
        status = False
        try:
            np.testing.assert_array_almost_equal(uv, uv_comp, decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            pass
        self.assertTrue(status)
    
    def test_get_rotation_operator(self):
        # Fixture
        u = np.random.rand(3)
        v = np.random.rand(3)
        
        # Test
        un = u / np.linalg.norm(u)
        vn = v / np.linalg.norm(v)

        R = util.get_rotation_operator(un, vn)
        unrot = R.dot(un)

        # Assert
        status = False
        try:
            np.testing.assert_array_almost_equal(vn, unrot, decimal=5)
            status = True
        except AssertionError as e:
            print e.message
            status = False
            pass
        self.assertTrue(status)
        
if __name__ == '__main__':
    unittest.main()