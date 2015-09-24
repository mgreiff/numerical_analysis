# -*- coding: utf-8 -*-

from numpy import array, linspace, dot
from numpy.linalg import solve
import matplotlib.pyplot as plt
import copy as copy
                
class Spline(object):
    """
    ARGS:
        dataX: A list containing float valued points.
        dataY: A list containing float valued points.
        intersect: A boolean True if the points are to be intersected by the
            B-splines, or False of the points describe a control polygon.
    """
    def __init__(self, dataX, dataY, intersect):
        self.curveX = []
        self.curveY = []
        self.minVal = 0.
        self.maxVal = 1.
        self.degree = 3 #Dysfunctional
        self.changePoints(dataX, dataY, intersect)

    def __call__(self, numberOfPoints):
        """
        Computes a B-spline curve for the currently stored controlX and -Y points.
        ARGS:
            numberOfPoints: The resolution of the caomputed B-spline (integer).
        """
        if not isinstance(numberOfPoints, int):
            raise ValueError('Unsupported "numberOfPoints" format %s. Must be of type int.' % str(type(numberOfPoints)))
        self.curveX = []
        self.curveY = []
        for element in linspace(.01, .99, numberOfPoints).tolist():
            i = (self.u > element).argmax() - 1
            self.curveX.append(self.d([None, None, None], element, i, self.controlX))
            self.curveY.append(self.d([None, None, None], element, i, self.controlY))

    def changePoints(self, dataX, dataY, intersect):
        """
        ARGS:
            dataX: A list containing float valued points.
            dataY: A list containing float valued points.
            intersect: A boolean True if the points are to be intersected by the
                B-splines, or False of the points describe a control polygon.
        """
        self.curveX = []
        self.curveY = []
        if not isinstance(dataX, list):
            raise ValueError('Unsupported dataX format %s. Must be of type list.' % str(type(dataX)))
        if not isinstance(dataY, list):
            raise ValueError('Unsupported dataY format %s. Must be of type list.' % str(type(dataY)))
        if not len(dataX) == len(dataY):
            raise Exception('Error. The data point arrays must be of equal length.')
        if not isinstance(intersect, bool):
            raise ValueError('Unsupported intersect option %s. Must be of type Boolean.' % str(type(intersect)))
        if intersect:
            self.pointsX = dataX
            self.pointsY = dataY
            self.u = array([self.minVal, self.minVal] +
                            linspace(self.minVal, self.maxVal, len(self.pointsX) - 2).tolist() +
                            [self.maxVal, self.maxVal])
            vander = self.getVander()           
            self.controlX = solve(vander, self.pointsX)
            self.controlY = solve(vander, self.pointsY)               
        else:
            self.controlX = dataX
            self.controlY = dataY
            self.u = array([self.minVal, self.minVal] +
                            linspace(self.minVal, self.maxVal, len(self.controlX) - 2).tolist() +
                            [self.maxVal, self.maxVal])
            vander = self.getVander() 
            self.pointsX = dot(vander, self.controlX)
            self.pointsY = dot(vander, self.controlY) 
    
    def getVander(self):
        """
        Returns a vandermonde matrix of basis dunctions evaluated at deBoor points
        points as calculated from the current u sequence.
        """
        xi = array([(self.u[i] + self.u[i + 1] + self.u[i + 2]) / 3 for i in range(len(self.u) - 2)])
        vander = array([[self.knot_sequence(self.u, i, self.degree)(xiVal) for i in range(len(xi))] for xiVal in xi])
        vander[len(xi) - 1][len(xi) - 1] = 1.
        return vander
        
    def plot(self, option):
        """
        Plots the results.
        ARGS:
            option: A string dscribing what the user wants to plot. Can either
                be set to of one of "basis", "control", "curve" or "XYpoints".
        """
        if not isinstance(option, str):
            raise ValueError('Unsupported "option" format %s. Must be of type str.' % str(type(option)))
        if option == 'basis':
            xx = linspace(0,1,100)
            for i in range(len(self.controlX)):
                yy = [self.knot_sequence(self.u, i, self.degree)(x) for x in xx]
                if i == len(self.controlX) - 1:
                    yy[99] = 1.0
                plt.plot(xx, yy, 'b')
        elif option == 'control':
            plt.plot(self.controlX, self.controlY, 'k-')
            plt.plot(self.controlX, self.controlY, 'ro', markersize=5)
        elif option == 'curve':
            if not self.curveX or not self.curveY:
                raise Exception('No curve has been generated. The __call__ method has to be run first.')
            plt.plot(self.curveX, self.curveY)
        elif option == 'XYpoints':
            plt.plot(self.pointsX, self.pointsY, 'gx', markersize=10)
            plt.plot(self.pointsX, self.pointsY, 'g-')
        else:
            raise TypeError('Unsupported option "%s". Must be of one of "basis", "control", "curve" or "XYpoints".' % option)

    def d(self, indices, uVal, i, dxy):
        if indices[0] != None:
            return dxy[indices[0]]
        elif not indices[2]:
            indicesA = [None, None, i + 1]
            indicesB = [None, None, i]
            return (self.alpha([i, i +1], uVal) * self.d(indicesB, uVal, i, dxy) + 
                    (1 - self.alpha([i, i + 1], uVal)) * self.d(indicesA, uVal, i, dxy))
        else:
            indexX  = 2 - [x for x in reversed(indices)].index(None)
            indicesA = copy.copy(indices)
            indicesA[indexX] = indices[indexX + 1]
            for ii in range(indexX + 1,2 + 1):
                indicesA[ii] = indicesA[ii] + 1
            indicesB = copy.copy(indices)
            indicesB[indexX] = indices[indexX + 1] - 1
            return (self.alpha(indicesA + indicesB, uVal) * self.d(indicesB, uVal, i, dxy) + 
                    (1 - self.alpha(indicesA + indicesB, uVal)) * self.d(indicesA, uVal, i, dxy))
        
    def alpha(self, arr, uVal):
        return (self.u[max(arr)] - uVal) / (self.u[max(arr)] - self.u[min([x for x in arr if x is not None])])

    def knot_sequence(self, p, i, k):
        if k == 0:
            if p[i - 1] == p[i]:
                return lambda u: 0
            else:
                return lambda u: 1 if p[i - 1] <= u and u < p[i] else 0
        else:
            if i + k >= len(p):
                return lambda u: ((u - p[i - 1]) / (p[i + k - 1] - p[i - 1]) *
                              self.knot_sequence(p, i, k - 1)(u))
            if p[i + k - 1] == p[i - 1]:
                return lambda u : ((p[i + k] - u)/(p[i + k] - p[i]) *
                              self.knot_sequence(p, i + 1, k - 1)(u))
            elif p[i + k] == p[i]:
                return lambda u: ((u - p[i - 1]) / (p[i + k - 1] - p[i - 1]) *
                              self.knot_sequence(p, i, k - 1)(u))
            else:
                return lambda u: ((u - p[i - 1]) / (p[i + k - 1] - p[i - 1]) *
                              self.knot_sequence(p, i, k - 1)(u) +
                              (p[i + k] - u)/(p[i + k] - p[i]) *
                              self.knot_sequence(p, i + 1, k - 1)(u))