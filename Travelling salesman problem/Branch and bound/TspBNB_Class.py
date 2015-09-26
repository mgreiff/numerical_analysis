# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

class TspBNB(object):
    
    def __init__(self, dataX, dataY):
        """
        Args:
            dataX (array): A 1xN array of x-points.
            dataY (array): A 1xN array of y-points.
        """
        self.xy = np.array([dataX, dataY])
        self.path = []
        self.fopt = None
        self.dist = self._distance_matrix()

    def __call__(self):
        print 'not yet implemented'

    def plot(self, points = 0, distance = 0, solution = 0):
        """
        Args:
            points (bool): Set totrue if the XY points are to be plotted
        """
        if solution:
            for ii in range(len(self.path) - 1):
                plt.gca().plot([self.xy[0][self.path[ii]], self.xy[0][self.path[ii + 1]]],
                               [self.xy[1][self.path[ii]], self.xy[1][self.path[ii + 1]]], 'k')
        if points:
            plt.gca().plot(self.xy[0], self.xy[1],'ro', markersize=5)         
            for i in range(len(self.xy[0])):
                plt.gca().annotate("#" + str(i), (self.xy[0][i],self.xy[1][i]))
                
        if distance:
            plt.gca().matshow(self.dist)
            

    def _distance_matrix(self):
        """
        Generates a distance matrix for the specified xy points.
        """
        nelm = len(self.xy[0])
        def dist(ii, jj):
            return (sqrt((self.xy[0][ii] - self.xy[0][jj]) ** 2 +
                   (self.xy[1][ii] - self.xy[1][jj]) ** 2))
        return np.array([np.array([dist(ii, jj) for jj in range(nelm)]) for ii in range(nelm)])

    def _travsalesman(self):
        nelm = len(self.xy[0])
        minmax = np.zeros([nelm,2]);
        for ii in range(nelm):
            minmax[ii][0]=min(np.concatenate((self.dist[:ii, ii], self.dist[ii + 1:, ii])));
            minmax[ii][1]=max(np.concatenate((self.dist[:ii, ii], self.dist[ii + 1:, ii])));
        startx = [0]
        bounds = self._boundy(startx, minmax)
        path, fopt = self._branchandbound(startx, minmax, bounds[1]);
        self.path = path
        self.fopt = fopt

    def _boundy(self, x, minmax):
        nelm = len(self.xy[0])
        bounds = np.array([0, 0])
        for ii in range(len(x) - 1):
            bounds = bounds + self.dist[x[ii]][x[ii + 1]]
        nonVisitedPoints = [ii for ii in range(nelm) if not ii in x]
        if nonVisitedPoints:
            for ii in nonVisitedPoints:
                bounds = bounds + minmax[ii][:]
        else:
            if len(x) <= nelm:
                bounds = bounds + minmax[x[0]][:]
        return bounds

    def _branchy(self, x):
        nelm = len(self.xy[0])
        if len(x) == nelm:
            possibleX = [x + [x[0]]]
        else:
            possibleX = []
            for ii in [ii for ii in range(nelm) if not ii in x]:
                possibleX.append(x + [ii])
        return possibleX

    def _branchandbound(self, x, minmax, fopt):
        bounds = self._boundy(x, minmax)
        if bounds[0] == bounds[1]:
            if bounds[1] < fopt:
                fopt = bounds[0]
        else:
            X = self._branchy(x)
            B = []
            for ii in range(len(X[:])):
                B.append(self._boundy(X[ii][:], minmax).tolist())
            for ii in [i[0] for i in sorted(enumerate(B), key=lambda x:x[1])]:
                if B[ii][0] < fopt:
                    xNew, foptNew = self._branchandbound(X[ii][:], minmax, fopt)
                    if foptNew < fopt:
                        fopt = foptNew
                        x = xNew
        return x, fopt