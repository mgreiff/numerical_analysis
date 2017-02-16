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

numberOfAnchors = 10
anchors = np.array([[1, 8, 8],
                    [0, 0, 0],
                    [2, 6, 3],
                    [2, 2, 3],
                    [6, 2, 3],
                    [5, 6, 7],
                    [0, 2, 4],
                    [2, 5, 3],
                    [2, 0, 3],
                    [3, 5, 0]])

# Here we assume perfect knowledge of the distances between the anchors
# this is our only data input to the calibration
Dsq = util.get_distance_matrix_squared(anchors)

# Define an anchor calibrator object
calibrator = AnchorCalibrator(numberOfAnchors)
calibrator.set_distance_matrix_squared(Dsq)

###############################
##### Specify constraints #####
###############################
# The index of the point which is to be the origin in the estimated coordinates
origin = 1
calibrator.set_origin(origin)
# A three nodes spanning a plane (here xy) used to determine the z-axis safely
plane = ((2,3,4),'xy')
calibrator.set_plane(plane)
# A two nodes residing on the same unit axis (here x)
axis = ((4,3),'x')
calibrator.set_axis(axis)

# Calibrate
calibrator()
estimatedPositions = calibrator.anchorPositions

# TODO add one more constraint to mirror correctly, two possible solutions
# always exist
estimatedPositions[:,1] *= -1.0;

###############################
### plot true configuration ###
###############################
fig = plt.figure(1)
ax = fig.add_subplot(121, projection = '3d')
sc = ax.scatter(anchors[:,0],
                anchors[:,1],
                anchors[:,2],
                edgecolor='red',
                linewidths=1)

for ii, txt in enumerate(np.arange(anchors.shape[0])):
    x2, y2, _ = proj3d.proj_transform(anchors[ii,0],anchors[ii,1],anchors[ii,2], ax.get_proj())
    plt.annotate(txt, xy = (x2, y2))

###############################
### plot generated solution ###
###############################
ax = fig.add_subplot(122, projection = '3d')

ax.scatter(estimatedPositions[:,0],
           estimatedPositions[:,1],
           estimatedPositions[:,2],
           edgecolor='green',
           linewidths=1)
for ii, txt in enumerate(np.arange(estimatedPositions.shape[0])):
    x2, y2, _ = proj3d.proj_transform(estimatedPositions[ii,0],
                                      estimatedPositions[ii,1],
                                      estimatedPositions[ii,2], ax.get_proj())
    plt.annotate(txt, xy = (x2, y2))

###############################
##### print the residual ######
###############################
residual = np.linalg.norm(estimatedPositions - anchors)
print "Residual accross all anchors: {0:.17f}".format(residual)
plt.show()