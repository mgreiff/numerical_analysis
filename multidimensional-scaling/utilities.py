# -*- coding: utf-8 -*-
from __future__ import division
 
import numpy as np
from math import sin, cos, acos
 
def multi_dimensional_scaling(Dsq):                                                                       
    numberOfPoints = Dsq.shape[0]                                                                    
    H = np.eye(numberOfPoints) - np.ones((numberOfPoints, numberOfPoints))/numberOfPoints                                                                                   
    B = -H.dot(Dsq).dot(H)/2                                                                            
    evals, evecs = np.linalg.eigh(B)                                                 
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]                   
    w =[0,1,2] # Three dimensions
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)
    return Y

def get_distance_matrix_squared(points):
    numberOfPoints = points.shape[0]
    Dsq = np.zeros((numberOfPoints,numberOfPoints))
    for ii in range(numberOfPoints):
        for jj in range(numberOfPoints):
            if ii != jj:
                Dsq[ii,jj] = np.linalg.norm(points[ii,:] - points[jj,:]) ** 2
    return Dsq

def get_skew_symmetric_operator(u):
    S = np.array([[0,     -u[2], u[1] ],
                  [u[2],  0,     -u[0]],
                  [-u[1], u[0],  0    ]])
    return S

def get_rotation_operator(a, b):
    v = np.cross(a, b)
    v = v/np.linalg.norm(v)
    S = get_skew_symmetric_operator(v)
    I = np.eye(3)
    theta = acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    R = I + sin(theta) * S + (1 - cos(theta)) * S.dot(S)    
    return R