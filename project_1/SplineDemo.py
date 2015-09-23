# -*- coding: utf-8 -*-

from SplineClass import Spline
import matplotlib.pyplot as plt

if __name__ == '__main__':

    numberOfPoints = 100 #The resolution of the spline curve
    if 1:
        # Generate spline by defining control points
        plt.figure(1)
        cX = [1, 2, 3, 5, 7, 8, 6]
        cY = [7, 4, 5, 1, 2, 4, 6]
        spl = Spline(cX, cY, False)
        spl(numberOfPoints)
        spl.plot('XYpoints')
        spl.plot('control')
        spl.plot('curve')
    if 1:
        # Generate spline and plot basis functions by defining x-, and y-points
        x = [0., 0.6, 1., 0.6, 0., 0., -0.6, -1., -0.6, 0.]
        y = [-1., 0., 0.9, 1.5, 1., 1., 1.5, 0.9, 0., -1.]
        spl = Spline(x, y, True)
        spl(numberOfPoints)
        plt.figure(2)
        spl.plot('XYpoints')
        #spl.plot('control')
        spl.plot('curve')
        plt.figure(3)
        spl.plot('basis')