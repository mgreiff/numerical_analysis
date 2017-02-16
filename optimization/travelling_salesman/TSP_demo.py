# -*- coding: utf-8 -*-
from TSP_classes import TspBNB, TspGA
import matplotlib.pyplot as plt
import random
import numpy as np

if __name__ == '__main__':
    print('\nThis example solves two TSP problems using the branch and\n'+
          'bound method and a genetic algorithm. The example is meant to be\n'+
          'run from the terminal, and the plots should be closed down to\n'
          'proceed when the program pauses\n')
          
    # Example implementation of the recursive TSP branch and bound algorithm
    # note that the problem may become computationally infeasible if N >~ 14
    N = 10                                                    # Number of elements to visit        
    x = np.array([random.randint(0,100) for ii in range(N)]) # Town x-coordinates
    y = np.array([random.randint(0,100) for ii in range(N)]) # Town y-coordinates
    tsp = TspBNB(x, y)                                       # TSP solver object
    tsp()                                                    # Call solver

    print('Visualizing results in problem 1 (BnB, N = %d randomized points)...' % (N))
    # Shows the distance matrix
    plt.figure()
    tsp.plot(distance = 1)
    plt.show()
    
    # Shows the towns and the computed optimal shortest route connecting them
    plt.figure()
    tsp.plot(solution = 1)
    plt.show()
    print('Complete!\n')
    
    # Example implementation of the recursive TSP branch and bound algorithm
    # note that the iteration limit or population size could be increased if
    # N >~ 100
    N = 50                                                   # Number of elements to visit 
    x = np.array([random.randint(0,100) for ii in range(N)]) # Town x-coordinates
    y = np.array([random.randint(0,100) for ii in range(N)]) # Town y-coordinates
    tsp = TspGA(x, y)                                        # TSP solver object
    tsp()

    print('Visualizing results in problem 2 (GA, N = %d randomized points)...' % (N))
    # Shows the distance matrix
    plt.figure()
    tsp.plot(distance = 1)
    plt.show()
    plt.figure()
    tsp.plot(solution = 1)
    plt.show()
    plt.figure()
    tsp.plot(solutionHistory = 1)
    plt.show()
    print('Complete!\n')