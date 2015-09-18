# -*- coding: utf-8 -*-

from numpy import array, linspace, zeros, append
from numpy.linalg import solve
import matplotlib.pyplot as plt
import copy as copy
import sys
                
class Spline(object):
    """
    ARGS:
        dataX: A list containing float valued points.
        dataY: --"--
        intersect: A boolean True if the points are to be intersected by the
            B-splines, or False of the points describe a control polygon.
    
    """
    def __init__(self, dataX, dataY, intersect):
        self.curveX = []
        self.curveY = []
        if not isinstance(intersect, bool):
            raise ValueError('Unsupported intersect option %s. Must be of type Boolean.' % str(type(intersect)))
        if intersect:
            print 'Not yet implemented'
        else:
            if not isinstance(dataX, list):
                raise ValueError('Unsupported dataX format %s. Must be of type list.' % str(type(dataX)))
            if not isinstance(dataY, list):
                raise ValueError('Unsupported dataY format %s. Must be of type list.' % str(type(dataY)))
            self.controlX = dataX
            self.controlY = dataY

    def __call__(self, numberOfPoints):
        if not isinstance(numberOfPoints, int):
            raise ValueError('Unsupported "numberOfPoints" format %s. Must be of type int.' % str(type(numberOfPoints)))
        self.curveX = []
        self.curveY = []
        uArr = array([0, 0] + linspace(0, 1, len(self.controlX) - 2).tolist() + [1, 1])
        for element in linspace(.01, .99, 50).tolist():
            i = (uArr > element).argmax() - 1
            self.curveX.append(self.d([None, None, None], element, i, self.controlX, uArr))
            self.curveY.append(self.d([None, None, None], element, i, self.controlY, uArr))

    def plot(self):
        plt.figure(1)
        plt.plot(self.controlX, self.controlY, color='k')
        plt.plot(self.controlX, self.controlY, 'ro', markersize=5)
        plt.plot(self.curveX, self.curveY)

    def d(self, indices, u, i, dxy, uArr):
        if indices[0] != None:
            return dxy[indices[0]]
        elif not indices[2]:
            indicesA = [None, None, i + 1]
            indicesB = [None, None, i]
            return (self.alpha([i, i +1], u, uArr) * self.d(indicesB, u, i, dxy, uArr) + 
                    (1 - self.alpha([i, i + 1], u, uArr)) * self.d(indicesA, u, i, dxy, uArr))
        else:
            indexX  = 2 - [x for x in reversed(indices)].index(None)
            indicesA = copy.copy(indices)
            indicesA[indexX] = indices[indexX + 1]
            for ii in range(indexX + 1,2 + 1):
                indicesA[ii] = indicesA[ii] + 1
            indicesB = copy.copy(indices)
            indicesB[indexX] = indices[indexX + 1] - 1
            return (self.alpha(indicesA + indicesB, u, uArr) * self.d(indicesB, u, i, dxy, uArr) + 
                    (1 - self.alpha(indicesA + indicesB, u, uArr)) * self.d(indicesA, u, i, dxy, uArr))
        
    def alpha(self, arr, u, uArr):
        return (uArr[max(arr)] - u) / (uArr[max(arr)] - uArr[min([x for x in arr if x is not None])])

if __name__ == '__main__':
    # Defines a set of control polygon points
    cX = [1, 2, 3, 5, 7, 8, 6]
    cY = [7, 4, 5, 1, 2, 4, 6]
    numberOfPoints = 50 #The resolution of the spline curve
    spl = Spline(cX, cY, False)
    if 0:
        spl.plot()
    if 1:
        spl(50)
        spl.plot()