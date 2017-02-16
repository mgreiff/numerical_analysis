# -*- coding: utf-8 -*-
import numpy as np
from numpy.matlib import repmat
import utilities as util

class AnchorCalibrator(object):

    def __init__(self, numberOfAnchors):
        self.numberOfAnchors = numberOfAnchors
        self.anchorPositions = None
        self.distanceMatrixSquared = None
        self.measurements = None
        self.origin = None
        self.axis = None
        self.plane = None

    def handle_measurements(self, measurements):
        """
        Takes a set of measurements from anchor i to anchor j with a
        distance d and computes the distance matrix
        """
        # Check std of eery measurment dimension
        # Remove outliers
        # Take mean distance and square the matrix
        print 'not yet written'
    
    def set_distance_matrix_squared(self, Dsq):
        self.distanceMatrixSquared = Dsq

    def __call__(self):

        if self.distanceMatrixSquared is None:
            raise Exception('Cannot compute solution as no distance matrix exists')

        # Computes the anchor positions
        self.anchorPositions = util.multi_dimensional_scaling(self.distanceMatrixSquared)

        # Rotates the solution if a plane has been given
        self.rotate_with_plane(self.plane)
            
        # Rotates the solution if an axis has been given
        self.rotate_with_axis(self.axis)

        # Change origin if an origin has been given
        self.change_origin(self.origin)

    def set_origin(self, origin):
        self.origin = origin

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
            temporaryOrigin = self.anchorPositions[ids[1],:]

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
            
            self.anchorPositions = self.anchorPositions - repmat(temporaryOrigin,self.numberOfAnchors,1)
            self.anchorPositions = R.dot(self.anchorPositions.T).T
            self.anchorPositions = self.anchorPositions + repmat(temporaryOrigin,self.numberOfAnchors,1)
            
    def rotate_with_axis(self, axis):
        if axis is not None:
            ids = self.axis[0]
            dimension = self.axis[1]
            u = self.anchorPositions[ids[0],:] - self.anchorPositions[ids[1],:]
            temporaryOrigin = self.anchorPositions[ids[1],:]

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
            self.anchorPositions = self.anchorPositions - repmat(temporaryOrigin,self.numberOfAnchors,1)
            self.anchorPositions = R.dot(self.anchorPositions.T).T
            self.anchorPositions = self.anchorPositions + repmat(temporaryOrigin,self.numberOfAnchors,1)

    def change_origin(self, origin):
        if origin is not None:
            shift = self.anchorPositions[origin,:]
            self.anchorPositions = self.anchorPositions - repmat(shift,self.numberOfAnchors,1)
