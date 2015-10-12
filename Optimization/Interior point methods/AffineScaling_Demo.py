# -*- coding: utf-8 -*-

from AffineScaling_Class import AffineScaling
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    x = np.array([0., 0.2])
    A = np.array([(2 * np.linspace(0,1,11)), np.ones(11)]).transpose()
    b = (np.linspace(0,1,11) ** 2 + 1)
    c = np.array([1., 1.])
    
    plt.figure(1)
    AffineScaling(x, A, b, c, plot = True, verbose = True)