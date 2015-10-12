# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class InteriorPointMethod(object):
    def __init__(self, x, A, b, c, plot = False, numericalLimit = 1e-5,
                 iterationLimit = 1e2, gamma = 0.5, verbose = False):

        self.A, self.x, self.b, self.c, self.gamma = A, x, b, c, gamma
        self.numericalLimit, self.iterationLimit = numericalLimit, iterationLimit
        self.plot, self.verbose, self.solutionHistory = plot, verbose, []
        self.__call__()
        
    def __call__(self, plot = False, verbose = False):
        self._sanity_check()
        if self.verbose: print 'Starting run... ',
        solution, numIter = self._iteration()
        if self.plot: self._plot([0, 1])
        if self.verbose: self._display(solution, numIter)
        return solution

    def _sanity_check(self):
        if np.shape(self.A)[0] != np.shape(self.b)[0] or np.shape(self.A)[1] != np.shape(self.c)[0]:
            raise Exception('The sized of the matrices do not match.')
        try:
            self.b = self.b.reshape(len(self.b),1)
            self.c = self.c.reshape(len(self.c),1)
            self.x = self.x.reshape(len(self.x),1)
        except:
            raise Exception('The arrays "b", "c" and "x" have to be one dimensional.')
        if self.gamma <= 0 or self.gamma >= 1:
            raise Exception('Gamma i outside of the interval (0,1).')

    def _iteration(self):
        print 'Needs to be re-implemented in specific interior point classes'

    def _plot(self, indices):
        # Plots function
        s0, s1 = self.solutionHistory[:, indices[0]], self.solutionHistory[:, indices[1]]
        resolution = 100
        x0min, x0max = min(s0) - (max(s0) - min(s0)) / 2., max(s0) + (max(s0) - min(s0)) / 2.
        x1min, x1max = min(s1) - (max(s1) - min(s1)) / 2., max(s1) + (max(s1) - min(s1)) / 2.
        x0 = np.linspace(x0min, x0max, resolution)
        x1 = np.linspace(x1min, x1max, resolution)
        Z = np.zeros((resolution, resolution))
        xHat = np.ones((len(self.x), 1))
        for ii in range(resolution):
            for jj in range(resolution):
                xHat[indices[0], 0] = x0[ii]
                xHat[indices[1], 0] = x1[jj]
                Z[ii, jj] = sum(np.dot(self.c.transpose(), xHat))
        X0, X1 = np.meshgrid(x0, x1)
        plt.colorbar(plt.contourf(X0, X1, Z, cmap=cm.gray), orientation='vertical')
        
        # Plots bounds
        for ii in range(np.shape(self.A)[0]):
            points = []
            xHat = (self.b[ii] - x0 * self.A[ii, indices[0]]) / self.A[ii, indices[1]]
            for jj in range(len(xHat)):
                if x1min < xHat[jj] and x1max > xHat[jj]:
                    points.append([x0[jj], xHat[jj]])
            if points:
                points = np.array(points).transpose()
                plt.plot(points[0], points[1], 'b')

        #Plots solution
        plt.plot(self.solutionHistory[:,indices[0]], self.solutionHistory[:,indices[1]], 'r')
        for ii in range(len(s0)): plt.plot(s0[ii], s1[ii], c='r', marker='.', markersize=10)
        plt.title(str(self)), plt.xlabel('$x_{%s}$' % indices[0]), plt.ylabel('$x_{%s}$' % indices[1])

    def _display(self, solution, numIter):
        print ('...run complete!\nSolution: %s\nNumber of iterations: %s\nNumerical limit: %s '
              % (solution.transpose()[0], numIter, self.numericalLimit))

class AffineScaling(InteriorPointMethod):

    def _iteration(self):
        k, x = 0, self.x
        self.solutionHistory = []
        while k < self.iterationLimit:
             v = self.b - np.dot(self.A, x)
             D = np.diag(v.reshape(1,11)[0] ** - 2)
             hx = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(self.A), D), self.A)), self.c)
             hv = -np.dot(self.A, hx)
             if hv.all() < 0:
                print 'Badly conditioned or unbounded problem'
             a = self.gamma * min(-v / hv)
             if self.plot: self.solutionHistory.append(x.reshape(1,2)[0])
             if np.linalg.norm(a * hx) < self.numericalLimit:
                 self.solutionHistory=  np.array(self.solutionHistory)
                 return x, k
             x = x + a * hx
             k += 1

    def __str__(self):
        return 'Affine scaling method'

class PrimalDual(InteriorPointMethod):
    # reference http://www.ieor.berkeley.edu/~ilan/ilans_pubs/affine_convex_1990.pdf
    # Complexity O(L\sqrt(n))
    def _iteration(self):
        k, x = 0, self.x
        z = self.c * self.x
        print z
        self.solutionHistory = []
        while k < self.iterationLimit:
            k += 1
        return x, k
        

    def __str__(self):
        return 'Primal-dual method'
x = np.array([0., 0.2])
A = np.array([(2 * np.linspace(0,1,11)), np.ones(11)]).transpose()
b = (np.linspace(0,1,11) ** 2 + 1)
c = np.array([1., 1.])

pd = PrimalDual(x, A, b, c, plot = False, verbose = True)