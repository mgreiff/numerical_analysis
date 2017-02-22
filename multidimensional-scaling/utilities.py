# -*- coding: utf-8 -*-
from __future__ import division
 
import numpy as np
from math import sin, cos, acos
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d

def multi_dimensional_scaling(Dsq):   
    dimension = 3                                                                    
    numberOfPoints = Dsq.shape[0]                                                                    
    J = np.eye(numberOfPoints) - np.ones((numberOfPoints, numberOfPoints))/numberOfPoints                                                                                   
    B = -J.dot(Dsq).dot(J)/2                                                                            
    evals, evecs = np.linalg.eigh(B)                                                 
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]                   
    Xsqrt  = np.diag(np.sqrt(evals[0:dimension]))
    E  = evecs[:,0:dimension]
    Y  = E.dot(Xsqrt)
    res = np.linalg.norm(evals[dimension:-1])
    return Y, res

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

def plot_anchors(points, figure, subplot, title, color):
    ax = figure.add_subplot(subplot, projection = '3d')
    ax.scatter(points[:,0], points[:,1], points[:,2], 'rx',edgecolor=color, linewidths=1)
    
    for ii, txt in enumerate(np.arange(points.shape[0])):
        x2, y2, _ = proj3d.proj_transform(points[ii,0],points[ii,1],points[ii,2], ax.get_proj())
        plt.annotate(txt, xy = (x2, y2))
    plt.title(title)