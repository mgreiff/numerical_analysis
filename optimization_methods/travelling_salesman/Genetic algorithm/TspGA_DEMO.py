# -*- coding: utf-8 -*-

from TspGA_Class import TspGA
import matplotlib.pyplot as plt
import numpy as np
import random as random

if __name__ == '__main__':
    # Example implementation of the genetic TSP algorithm
    if 1:
        N = 50
        x = np.array([random.randint(0,N) for ii in range(N)])
        y = np.array([random.randint(0,N) for ii in range(N)])
        tsp = TspGA(x, y)
        plt.figure(1)
        tsp.plot(distance = 1)
        plt.figure(2)
        tsp._genetic()
        tsp.plot(solution = 1, points = 1)
        plt.figure(3)
        tsp.plot(solutionHistory = 1)
