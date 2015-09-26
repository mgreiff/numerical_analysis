# -*- coding: utf-8 -*-

from TspGA_Class import TspGA
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # Example implementation of the genetic TSP algorithm
    x = np.array([0., 0.6, 1., 0.6, 0., -0.6, -1., -0.6,])
    y = np.array([-1., 0., 0.9, 1.5, 1., 1.5, 0.9, 0.,])
    tsp = TspGA(x, y)
    if 1:
        plt.figure(1)
        tsp.plot(distance = 1)
    if 0:
        plt.figure(2)
        tsp.plot(solution = 1, points = 1)
    if 0:
        plt.figure(3)
        x = np.array([1, 2, 3, 5, 7, 8, 6])
        y = np.array([7, 4, 5, 1, 2, 4, 6])
        tsp = TspGA(x, y)
        tsp.plot(solution = 1, points = 1)
