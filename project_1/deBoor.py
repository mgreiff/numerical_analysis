# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 19:50:45 2015

@author: Mgreiff
"""

from numpy import array, linspace, zeros
from numpy.linalg import solve
import matplotlib.pyplot as plt

class Spline(object):
    def __init__(self, points, k = 3, i = 1):
        self.points = points
        self.k = k
        self.i = i

    def __call__(self, *args):
        if len(args) == 1:
            self.k = args[0]
        elif len(args) > 1:
            self.k, self.i = args[0], args[1]
        return self.knot_sequence(self.points, self.i, self.k)

    def plot(self):
        return 'Not yet implemented'
        
    def knot_sequence(self, p, i, k):
        """
        Creates the basis functions of the interpolation recursively, calling
        itself until k == 0, where the function N^0_i is defined as described
        on page 1.3 in unit 1.
        ARGS:
            p: An array of floats of with u-values [u_0, u_1, ..., u_K].
            i: The integer index i of an element in p (i.e., u_i) to which the
                basis function corresponds (this is the first element in p at
                the function is not equal to zero).
            k: The integer describing the degree of the created basis function.
        """
        if k == 0:
            if p[i - 1] == p[i]:
                return lambda u: 0
            else:
                return lambda u: 1 if p[i - 1] <= u and u < p[i] else 0
        else:
            return lambda u: ((u - p[i - 1]) / (p[i + k - 1] - p[i - 1]) *
                              self.knot_sequence(p, i, k - 1)(u) +
                              (p[i + k] - u)/(p[i + k] - p[i]) *
                              self.knot_sequence(p, i + 1, k - 1)(u))

if __name__ == '__main__':
    if 1:
        # Creates 7 points (marked red) in the interval [0,1], plots them and the
        # basis unctions u_1 (blue) and u_2 green.
        points = linspace(0, 1, 7).tolist()
        spl = Spline(points)                          # Creates a spline instance
        xx = linspace(min(points), max(points), 100)  # Defines the x grid
        y1 = [spl(3, 1)(x) for x in xx]               # Calculates the u_1 basis f
        y2 = [spl(3, 2)(x) for x in xx]               # --"--
        plt.figure(1)
        plt.plot(xx, y1, 'b')
        plt.plot(xx, y2, 'g')
        plt.plot(points, zeros(len(points)), 'ro', markersize=5)

    if 0:
        x = [1, 2, 3, 5, 7, 8, 7]
        y = [7, 4, 5, 1, 2, 4, 6]
        u = linspace(0, 7, len(x) + 2).tolist()
        spl = Spline(u)
        xi = array([(u[i] + u[i + 1] + u[i + 2]) / 3 for i in range(len(u) - 2)])
        #Obs! indexeringen funkar inte riktigt här
        #Vandermonde = array([[spl(3,i)(xiVal) for i in range(len(x))] for xiVal in xi])
        plt.figure(2)
        plt.plot(x, y, color='k')
        plt.plot(x, y, 'ro', markersize=10)
    
    