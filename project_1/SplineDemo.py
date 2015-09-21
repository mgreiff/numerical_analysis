# -*- coding: utf-8 -*-

from SplineClass import Spline
import matplotlib.pyplot as plt

if __name__ == '__main__':

    numberOfPoints = 50 #The resolution of the spline curve
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
        x = [1, 2, 3, 5, 7, 8, 7]
        y = [7, 4, 5, 1, 2, 4, 6]
        spl = Spline(x, y, True)
        spl(numberOfPoints)
        plt.figure(2)
        spl.plot('XYpoints')
        spl.plot('control')
        spl.plot('curve')
        plt.figure(3)
        spl.plot('basis')
    if 1:
        uval = 0.4
        i = (spl.u > uval).argmax() - 1
        print spl.d([None, None, None], 0.4, i, spl.controlX, spl.u)
        print sum([spl.knot_sequence(spl.u, i, 3)(uval) * spl.controlX[i] for i in range(len(spl.controlX))])