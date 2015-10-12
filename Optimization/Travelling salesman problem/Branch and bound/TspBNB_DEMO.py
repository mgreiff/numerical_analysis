# -*- coding: utf-8 -*-

from TspBNB_Class import TspBNB
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # Example implementation of the recursive TSP algorithm
    x = np.array([0., 0.6, 1., 0.6, 0., -0.6, -1., -0.6,])
    y = np.array([-1., 0., 0.9, 1.5, 1., 1.5, 0.9, 0.,])
    tsp = TspBNB(x, y)
    if 1:
        plt.figure(1)
        tsp.plot(distance = 1)
    if 1:
        plt.figure(2)
        tsp._travsalesman()
        tsp.plot(solution = 1, points = 1)
    if 1:
        plt.figure(3)
        x = np.array([1, 2, 3, 5, 7, 8, 6])
        y = np.array([7, 4, 5, 1, 2, 4, 6])
        tsp = TspBNB(x, y)
        tsp._travsalesman()
        tsp.plot(solution = 1, points = 1)
