# -*- coding: utf-8 -*-
"""
This example demonstrates how solutions to the TOA calibration
problem may be found using multi-dimensional scaling. The example
sets 10 anchors in space and the recreates their position based
on their mutual distances and a set of constraints
"""
import matplotlib.pylab as plt
import numpy as np
import utilities as util
from classes import AnchorCalibrator
from mpl_toolkits.mplot3d import Axes3D, proj3d
import sys

# True data
numberOfAnchors = 20
anchors = np.array([[1, 8, 8],
                    [0, 0, 0],
                    [2, 6, 3],
                    [2, 2, 3],
                    [6, 2, 3],
                    [5, 6, 7],
                    [0, 2, 4],
                    [2, 5, 3],
                    [2, 4, 3],
                    [3, 5, 6],
                    [2, 9, 9],
                    [1, 1, 1],
                    [3, 7, 4],
                    [3, 3, 4],
                    [7, 3, 4],
                    [6, 7, 8],
                    [1, 3, 5],
                    [3, 6, 4],
                    [3, 5, 4],
                    [4, 6, 5]])


# Define an anchor calibrator object
calibrator = AnchorCalibrator(numberOfAnchors)

if len(sys.argv) == 1:
    print 'Provide a string argument, either "true" or "noisy".'
else:
    if sys.argv[1] == 'true':
        # Here we assume perfect knowledge of the distances between the anchors
        Dsq = util.get_distance_matrix_squared(anchors)
        calibrator.set_distance_matrix_squared(Dsq)
    elif sys.argv[1] == 'noisy':
        # Include measurements as the come
        for ii in range(10000):
            i = np.random.randint(20)
            j = np.random.randint(20)
            if i != j:
                noise = np.random.normal(0.01,0.01)
                distance = np.linalg.norm(anchors[i,:] - anchors[j,:]) + noise
                dataPoint = ((i, j), distance)
                # include new data
                calibrator.include_data_point(dataPoint)
        # Process data
        calibrator.process_data()

    # The index of the point which is to be the origin in the estimated coordinates
    origin = 1
    calibrator.set_origin(1)
    
    # A three nodes spanning a plane (here xy) used to determine the z-axis safely
    plane = ((2,3,4),'xy') # from 2 to 3
    calibrator.set_plane(plane)
    
    # A two nodes residing on the same unit axis (here x)
    axis = ((3,2),'y') # from index 3 to index 2
    calibrator.set_axis(axis)
    
    # mirror
    normal = ((3,4),'x') # from index 3 to index 4
    calibrator.set_mirror(normal)
    
    # Calibrate
    calibrator()
    
    # Extract estimate
    estimate = calibrator.anchorPositions
    
    # TODO add one more constraint to mirror correctly
    #estimate[:,0] *= -1.0;
    
    ###############################################################################
    ### plot true and estimated anchors, compute residual                       ###
    ###############################################################################
    fig = plt.figure(1)
    util.plot_anchors(anchors,
                      fig, 121, 'True anchor positions', 'red')
    util.plot_anchors(estimate,
                      fig, 122, 'Estimated anchor positions', 'green')
    
    print "Residual accross all anchors: {0:.17f}".format(np.linalg.norm(estimate - anchors))
    plt.show()