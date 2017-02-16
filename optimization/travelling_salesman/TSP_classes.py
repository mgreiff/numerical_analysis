# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import random
import copy as copy
import sys

class TspBNB(object):
    """
    The branch and bound TSP object, enabling the setting of positions in the
    constructor, a method for plotting
    """
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

    def plot(self, solution = 0, distance = 0):
        """
        Args:
            solution (bool): Set to true if the XY points and solution are to
                be plotted
            distance (bool): Set to true if the distance matrix is to be
                plotted
        """
        if solution:
            for ii in range(len(self.path) - 1):
                plt.gca().plot([self.xy[0][self.path[ii]], self.xy[0][self.path[ii + 1]]],
                               [self.xy[1][self.path[ii]], self.xy[1][self.path[ii + 1]]], 'k')
            plt.gca().plot(self.xy[0], self.xy[1],'ro', markersize=5) 
            for i in range(len(self.xy[0])):
                plt.gca().annotate("#" + str(i), (self.xy[0][i],self.xy[1][i]))
            plt.title(('Optimal solution to the TSP (black) with connected towns,\n'+
                       '(red) using the branch and bound method'))
                        
        if distance:
            ax = plt.gca()
            cax = ax.matshow(self.dist, interpolation='nearest')
            plt.title('Distance matrix $\in\mathbb{R}^{NxN}$ from town i to town j')
            plt.colorbar(cax)
            

    def _distance_matrix(self):
        """
        Generates and returns a distance matrix for the specified xy points.
        """
        nelm = len(self.xy[0])
        def dist(ii, jj):
            return (sqrt((self.xy[0][ii] - self.xy[0][jj]) ** 2 +
                   (self.xy[1][ii] - self.xy[1][jj]) ** 2))
        return np.array([np.array([dist(ii, jj) for jj in range(nelm)]) for ii in range(nelm)])

    def __call__(self):
        """
        Solves the traveling salesman problem using information set in the
        constructor.
        """
        sys.stdout.write('Starting BnB TSP solver... '),
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
        sys.stdout.write('Complete!\n')
        sys.stdout.flush()

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

class TspGA(object):
    
    def __init__(self, dataX, dataY, ):
        """
        Args:
            dataX (array): A 1xN array of x-points.
            dataY (array): A 1xN array of y-points.
        """
        self.xy = np.array([dataX, dataY])
        self.populationSize = 100
        self.iterationLim = 1000
        self.solutionHistory = []
        self.path = []
        self.dist = self._distance_matrix()

    def plot(self, points = 0, distance = 0, solution = 0, solutionHistory = 0):
        """
        Args:
            points (bool): Set to true if the XY points are to be plotted
        """
        if solution:
            for ii in range(len(self.path) - 1):
                plt.gca().plot([self.xy[0][self.path[ii]], self.xy[0][self.path[ii + 1]]],
                               [self.xy[1][self.path[ii]], self.xy[1][self.path[ii + 1]]], 'k')
            plt.gca().plot(self.xy[0], self.xy[1],'ro', markersize=5)
            plt.xlabel('x'), plt.ylabel('y')
            plt.title('Approximate solution to the TSP (black) with connected towns,\n'+
                       '(red) using the branch and bound method')
            for i in range(len(self.xy[0])):
                plt.gca().annotate("#" + str(i), (self.xy[0][i],self.xy[1][i]))
        if distance:
            plt.gca().matshow(self.dist)
            plt.title('Distance matrix $\in\mathbb{R}^{NxN}$ from town i to town j')
        if solutionHistory:
            plt.gca().plot([ii for ii in range(self.iterationLim)], self.solutionHistory)
            plt.title('Path distance (cost) as a function of the iterationnumber')
            plt.ylabel('Cost'), plt.xlabel('Number of iterations')
        plt.show()

    def _distance_matrix(self):
        """
        Generates a distance matrix for the specified xy points.
        """
        def dist(ii, jj):
            """
            Calculates a distance between two points at indices ii and jj in
            the xy data matrix.
            ARGS:
                ii, jj (int): Indices
            """
            return (sqrt((self.xy[0][ii] - self.xy[0][jj]) ** 2 + (self.xy[1][ii] - self.xy[1][jj]) ** 2))
        return np.array([np.array([dist(ii, jj) for jj in range(len(self.xy[0]))]) for ii in range(len(self.xy[0]))])
    
    def __call__(self):
        # Generates a population
        print('Starting GA solver...')
        n = len(self.xy[0])
        population = np.array([[ii for ii in range(n)] for jj in range(self.populationSize)])
        for ii in range(1, self.populationSize):
            random.shuffle(population[ii])

        for iteration in range(self.iterationLim):
            self.print_progress(iteration, self.iterationLim, 'Status:', 'complete.')
            # Computes the total distance for each population member
            populationDist = np.array([sum([self.dist[population[jj][ii - 1], population[jj][ii]] for ii in range(n)]) for jj in range(self.populationSize)])
            randomizedIndices = [ii for ii in range(self.populationSize)]
            random.shuffle(randomizedIndices)
            
            path, minDistance = self._bestSolution(population)
            self.path = path + [path[0]]
            self.solutionHistory.append([minDistance])
        
            newPopulation = []
            for ii in range(self.populationSize // 4):
                selectedPopulations = population[randomizedIndices[4 * ii : 4 * (ii + 1)]]
                selectedDistances = populationDist[randomizedIndices[4 * ii : 4 * (ii + 1)]]
                index = np.where(selectedDistances == selectedDistances.min())[0][0]
                bestRoute = selectedPopulations[index]
                breakPoints = [random.randint(0, n - 1), random.randint(0, n - 1)]
                breakPoints.sort(key=int)

                offspring = copy.copy(bestRoute)
                    
                flip = copy.copy(bestRoute).tolist()    #Flip
                flipSection = flip[breakPoints[0]:breakPoints[1]]
                flipSection.reverse()
                flip = flip[:breakPoints[0]] + flipSection + flip[breakPoints[1]:]
                offspring = np.append([offspring], [flip], axis=0)
                         
                swap = copy.copy(bestRoute)             #Swap
                swap[breakPoints[0]], swap[breakPoints[1]] = swap[breakPoints[1]], swap[breakPoints[0]]
                offspring = np.append(offspring, [swap], axis=0)
            
                slide = copy.copy(bestRoute).tolist()   #Slide
                poppedElement = slide.pop(breakPoints[1])
                slide.insert(breakPoints[0], poppedElement)
                offspring = np.append(offspring, [slide], axis=0)
            
                if newPopulation == []:
                    newPopulation = offspring
                else:
                    newPopulation = np.append(newPopulation, offspring, axis=0)
            population = newPopulation
    
    def _bestSolution(self, population):
        populationDist = np.array([sum([self.dist[population[jj][ii - 1], population[jj][ii]] for ii in range(len(population[0]))]) for jj in range(self.populationSize)])
        index = np.where(populationDist == populationDist.min())
        return population[index[0][0]].tolist(), populationDist.min()

    def print_progress(self, iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 30):
        """
        Prints progress bar in terminal window
        Args:
           iteration - positive integer. The current iteration.
           total - positive integer > iteration. The total number of iterations before completion.
           prefix - string. Empty by default, specifies text before the progress bar.
           suffix - string. Empty by default, specifies text after the progress bar.
           decimals - positive integer. number of decimals in the percentage calculation.
           barLength - positive, non-zero integer. Set to 30 # by default.
        """
        filledLength = int(round(barLength * iteration / float(total)))
        percents = round(100.00 * (iteration / float(total)), decimals)
        bar = '#' * filledLength + '-' * (barLength - filledLength)
        sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
        sys.stdout.flush()
        if iteration == total-1:
            sys.stdout.write('Complete!' + ' ' * (barLength + 20) + '\n')