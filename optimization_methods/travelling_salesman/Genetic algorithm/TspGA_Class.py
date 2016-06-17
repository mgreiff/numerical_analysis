# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import random
import copy as copy

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

    def __call__(self):
        print 'not yet implemented'

    def plot(self, points = 0, distance = 0, solution = 0, solutionHistory = 0):
        """
        Args:
            points (bool): Set to true if the XY points are to be plotted
        """
        if solution:
            for ii in range(len(self.path) - 1):
                plt.gca().plot([self.xy[0][self.path[ii]], self.xy[0][self.path[ii + 1]]],
                               [self.xy[1][self.path[ii]], self.xy[1][self.path[ii + 1]]], 'k')
            plt.xlabel('x'), plt.ylabel('y'), plt.title('Solution')
        if points:
            plt.gca().plot(self.xy[0], self.xy[1],'ro', markersize=5)
            plt.xlabel('x'), plt.ylabel('y'), plt.title('Input points')
            for i in range(len(self.xy[0])):
                plt.gca().annotate("#" + str(i), (self.xy[0][i],self.xy[1][i]))
        if distance:
            plt.gca().matshow(self.dist)
            plt.title('Distance matrix')
        if solutionHistory:
            plt.gca().plot([ii for ii in range(self.iterationLim)], self.solutionHistory)
            plt.title('Solution history')
            plt.ylabel('Cost'), plt.xlabel('Number of iterations')
        if solution and points:
            plt.title('Input points and solution')
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
    
    def _genetic(self):
        # Generates a population
        n = len(self.xy[0])
        population = np.array([[ii for ii in range(n)] for jj in range(self.populationSize)])
        for ii in range(1, self.populationSize):
            random.shuffle(population[ii])

        for iteration in range(self.iterationLim):
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
                    
                flip = copy.copy(bestRoute).tolist() #Flip
                flipSection = flip[breakPoints[0]:breakPoints[1]]
                flipSection.reverse()
                flip = flip[:breakPoints[0]] + flipSection + flip[breakPoints[1]:]
                offspring = np.append([offspring], [flip], axis=0)
                         
                swap = copy.copy(bestRoute) #Swap
                swap[breakPoints[0]], swap[breakPoints[1]] = swap[breakPoints[1]], swap[breakPoints[0]]
                offspring = np.append(offspring, [swap], axis=0)
            
                slide = copy.copy(bestRoute).tolist() #Slide
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

