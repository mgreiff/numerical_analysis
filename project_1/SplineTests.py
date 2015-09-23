# -*- coding: utf-8 -*-

from SplineClass import Spline
import unittest
import random

class TestStringMethods(unittest.TestCase):
    def setUp(self):
        length = random.randint(5,10)
        self.x = [random.randint(0,10) for i in range(length)]
        self.y = [random.randint(0,10) for i in range(length)]
        self.spl = Spline(self.x, self.y, True)
        self.numberOfPoints = random.randint(30, 50)
        self.spl(self.numberOfPoints)
    
    def testSumOfBaseFunctions(self):
        uval = random.random()
        a = sum([self.spl.knot_sequence(self.spl.u, i, 3)(uval) for i in range(len(self.spl.controlX))])
        self.assertAlmostEqual(a, 1.)
    
    def testMultiply(self):
        uval = random.random()
        index = (self.spl.u > uval).argmax() - 1
        a = self.spl.d([None, None, None], uval, index, self.spl.controlX)
        b = sum([self.spl.knot_sequence(self.spl.u, i, 3)(uval) * self.spl.controlX[i] for i in range(len(self.spl.controlX))])
        self.assertAlmostEqual(a, b)

if __name__ == '__main__':
    unittest.main()