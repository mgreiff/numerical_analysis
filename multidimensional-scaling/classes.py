# -*- coding: utf-8 -*-
import numpy as np
from numpy.matlib import repmat
import utilities as util

class AnchorCalibrator(object):

    def __init__(self, numberOfAnchors):
        # Public attributes
        self.numberOfAnchors = numberOfAnchors
        self.anchorPositions = None

        # Private attributes
        self.distanceMatrixSquared = None
        self.data = None
        self.means = None
        self.deviations = None
        self.reset_data()

        # Constraints
        self.origin = None
        self.axis = None
        self.normal = None
        self.plane = None

    def include_data_series(self, measurements):
        """
        Takes a set of measurements from anchor i to anchor j with a
        distance d and computes the distance matrix. Each measurement is
        represented as a tuple, ((i,j),distance).
        """
        for measurement in measurements:
            self.include_data_point(measurement)

    def include_data_point(self, measurement):
        """
        Takes a set of measurements from anchor i to anchor j with a
        distance d and computes the distance matrix. Each measurement is
        represented as a tuple, ((i,j),distance).
        """
        i = measurement[0][0]
        j = measurement[0][1]
        distance = measurement[1]
        if i > j:
            self.data[i][j].append(distance)
        else:
            self.data[j][i].append(distance)

    def reset_data(self):
        self.means = np.zeros((self.numberOfAnchors,self.numberOfAnchors))
        self.deviations = np.zeros((self.numberOfAnchors,self.numberOfAnchors))
        self.distanceMatrixSquared = np.zeros((self.numberOfAnchors,self.numberOfAnchors))
        self.data = [[ [] for ii in range(self.numberOfAnchors)] for jj in range(self.numberOfAnchors)]

    def process_data(self):
        for ii in range(self.numberOfAnchors) :
            for jj in range(self.numberOfAnchors):
                if ii > jj:
                    # Compute mean and standard deviation
                    if not self.data[ii][jj]:
                        print 'here'
                        self.means[ii,jj] = np.nan
                        self.deviations[ii,jj] = np.nan
                    else:
                        measurements_ij = np.array(self.data[ii][jj])
                        self.means[ii,jj] = np.mean(measurements_ij)
                        self.deviations[ii,jj] = np.std(measurements_ij)
        meanSquared = self.means ** 2
        # Do some outlier rejection
        self.distanceMatrixSquared = meanSquared + meanSquared.T
        
    def __call__(self):
        if self.distanceMatrixSquared is None:
            raise Exception('Cannot compute solution as no distance matrix exists')

        # Computes the anchor positions
        self.anchorPositions, residual = util.multi_dimensional_scaling(self.distanceMatrixSquared)

        # Rotates the solution if a plane has been given
        self.rotate_with_plane(self.plane)
            
        # Rotates the solution if an axis has been given
        self.rotate_with_axis(self.axis)

        # Mirrors the system
        self.mirror_with_axis(self.normal)

        # Change origin if an origin has been given
        self.change_origin(self.origin)

    def set_distance_matrix_squared(self, Dsq):
        self.distanceMatrixSquared = Dsq

    def set_origin(self, origin):
        self.origin = origin

    def set_mirror(self, normal):
        self.normal = normal

    def set_axis(self, axis):
        self.axis = axis

    def set_plane(self, plane):
        self.plane = plane
    
    def rotate_with_plane(self, plane):
        if plane is not None:
            ids = self.plane[0]
            dimension = self.plane[1]
            u = self.anchorPositions[ids[0],:] - self.anchorPositions[ids[1],:]
            v = self.anchorPositions[ids[2],:] - self.anchorPositions[ids[1],:]

            # Current orientation
            uv = np.cross(u, v)
            current = uv / np.linalg.norm(uv)

            #Intended orientation
            if not 'x' in dimension:
                intended = np.array([1,0,0])
            if not 'y' in dimension:
                intended = np.array([0,1,0])
            if not 'z' in dimension:
                intended = np.array([0,0,1])
            R = util.get_rotation_operator(current, intended)

            self.anchorPositions = R.dot(self.anchorPositions.T).T

            
    def rotate_with_axis(self, axis):
        if axis is not None:
            ids = self.axis[0]
            dimension = self.axis[1]
            u = self.anchorPositions[ids[1],:] - self.anchorPositions[ids[0],:]

            # Current orientation
            current = u / np.linalg.norm(u)
            
            #Intended orientation
            if 'x' in dimension:
                intended = np.array([1,0,0])
            if 'y' in dimension:
                intended = np.array([0,1,0])
            if 'z' in dimension:
                intended = np.array([0,0,1])

            R = util.get_rotation_operator(current, intended)
            self.anchorPositions = R.dot(self.anchorPositions.T).T
            
            u = self.anchorPositions[ids[0],:] - self.anchorPositions[ids[1],:]

    def mirror_with_axis(self, normal):
        if normal is not None:
            ids = self.normal[0]
            dimension = self.normal[1]
            ntilde = self.anchorPositions[ids[1],:] - self.anchorPositions[ids[0],:]

            # Current orientation
            ntilde = ntilde / np.linalg.norm(ntilde)
            
            #Intended orientation
            if 'x' in dimension:
                n = np.array([1,0,0])
                col = 0
            if 'y' in dimension:
                n = np.array([0,1,0])
                col = 1
            if 'z' in dimension:
                n = np.array([0,0,1])
                col = 2
            product = n.dot(ntilde)

            if product < 0.0:
                self.anchorPositions[:,col] *= -1.0
    
    def change_origin(self, origin):
        if origin is not None:
            shift = self.anchorPositions[origin,:]
            self.anchorPositions = self.anchorPositions - repmat(shift,self.numberOfAnchors,1)
