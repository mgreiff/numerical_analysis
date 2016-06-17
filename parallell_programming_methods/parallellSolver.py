# -*- coding: utf-8 -*-
from mpi4py import MPI
import sys
import matplotlib.pyplot as plt
import json
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

class Room(object):
    def __init__(self, comm, deltax):
        """
        Initalizes grid, sets the boundaries for the rooms in the apartment
        and interconnects the rooms to each other.
        ArgsARGS:
            comm [MPI object]: communcicator between processes
            deltax [float]: distance between adjacent nodes
        """
        self.comm = comm
        self.rank = comm.Get_rank()
        self.iterMax = 10;
        self.neumannBC = []

        print 'Initializing process %s ...' % self.rank
        self.deltax = deltax
        self.M = int(round(1/deltax))
        self.N = int(round(1/deltax))

        if self.rank == 0:
            bc = [15., 15., 40., None]
            self.indices, self.boundaries = self.initializeGrid(bc)
            self.neighbours = {'rooms' : [1], 'indices' : (np.arange(20)*20 - 1).tolist()[2:]}
        elif self.rank == 1:
            self.M = 2*self.M - 1
            bc = [40., 5., 15., 15.]
            self.indices, self.boundaries = self.initializeGrid(bc)
            self.neighbours = {'rooms' : [0, 2],
                                'indices' : [(400 + np.arange(20)*20).tolist()[:-2], (np.arange(20)*20 - 1).tolist()[2:]]}
        elif self.rank == 2:  
            bc = [15., 15., None, 40.]
            self.indices, self.boundaries = self.initializeGrid(bc)
            self.neighbours = {'rooms' : [1], 'indices' : (np.arange(20)*20).tolist()[1:-1]}
        else:
            raise ValueError('Faulty value: ' + str(self.rank))

    def tempDistribution(self):
        """
        Comutes the temperature distribution for the apartment after
        self.iterMax iterations.
        Return:
            solution [numpy array (M x N)]: Returns a matrix containing
            temperatures at the nodes.
        """
        for iteration in range(self.iterMax):
            for destination in self.neighbours['rooms']:
                if not (iteration == 0 and self.rank == 1):
                    data = self.comm.recv(source=destination)
                    if self.rank == 1:
                        ind = self.neighbours['rooms'].index(destination)
                        gridIndex = self.neighbours['indices'][ind]
                        for ii in range(len(gridIndex)):
                            self.boundaries[self.indices.index(gridIndex[ii])] = data[ii]
                    else:
                        self.neumannBC = data

            newSolution = self.solveSystem()
            if iteration == 0:
                solution = newSolution
            else:
                solution = 0.8*newSolution + 0.2*solution
                
            if self.rank != 1:
                sendData = solution[self.neighbours['indices']]
                self.comm.send(sendData, dest=self.neighbours['rooms'][0])
            else:
                sR = solution[[ii + 1 for ii in self.neighbours['indices'][0]]]
                sendData1 = (sR - solution[self.neighbours['indices'][0]]) / self.deltax
                self.comm.send(sendData1, dest=0)

                sL = solution[[ii - 1 for ii in self.neighbours['indices'][1]]]
                sendData2 = (sL - solution[self.neighbours['indices'][1]])  / self.deltax
                self.comm.send(sendData2, dest=self.neighbours['rooms'][1])
        print '...iteration in rank %s complete!' % self.rank
        return np.reshape(solution, [self.M, self.N])

    def solveSystem(self):
        """
        Solves the 2D-laplace equation on a rectangular grid given the
        information stored in the Room object.
        Return:
            solution [np.array (M x N)]: The temperatures at all gridpoints,
                represented as a matrix (array of arrays).
        """
        M = self.M
        N = self.N
        deltax = self.deltax
        indices= self.indices
        boundaries = self.boundaries
        neumannindices = []

        if self.rank != 1:
            neumannindices = self.neighbours['indices']
            neumann = self.neumannBC

        unknowns = list(set(np.arange(M * N).tolist()) - set(self.indices))
        nelm = len(unknowns)
        A = np.zeros((nelm, nelm))
        b = np.zeros((nelm, 1))

        for val, ii in zip(unknowns, range(len(unknowns))):
            if val in neumannindices:
                A[ii,ii] = -3.
                b[ii] = -neumann[neumannindices.index(val)] * (deltax)
                if self.rank == 0:
                    surroundingValues = [val - 1, val - N, val + N] #[left, up, down]
                if self.rank == 2:
                    surroundingValues = [val + 1, val - N, val + N] #[right, up, down]
            else:
                A[ii,ii] = -4.
                surroundingValues = [val + 1, val - 1, val - N, val + N] #[right, left, up, down]
  
            for element in surroundingValues:
                if element in indices:
                    b[ii] -= boundaries[indices.index(element)]
                else:
                    A[ii, unknowns.index(element)] = 1.

        u = np.linalg.solve(A, b)
        solution = np.zeros(M*N)
        for ii in range(M):
            for jj in range(N):
                ind = ii * N + jj
                if ind in self.indices:
                    solution[ind] = self.boundaries[self.indices.index(ind)]
                else:
                    solution[ind] = u[unknowns.index(ind)]
        return solution

    def initializeGrid(self, bc):
        """
        Initializes an M x N grid with boundary conditions bc. The boundary
        conditions are set to bc[0], bc[0], bc[0], bc[3] along the boundaries
        on the North, South, West and East side of the grid. If an element in
        bc is set to None, no boundary condition is set along the corresponding
        boundary.
        ARGS:
            bc [list]: The boundary conditions along the sides of the
                rectangular grid.
        Return:
            indices [list]: A list containing indices of elements in the grid.
            boundaries [list]: A list with temperatures of gridpoints on 
                defined boundaries.
        """
        indices = []
        boundaries = []
        M = self.M
        N = self.N

        if bc[0] != None:
            for ii in range(N):
                boundaries.append(bc[0])
                indices.append(ii)
        if bc[1] != None:
            for ii in range(M*N - N, M*N):
                boundaries.append(bc[1])
                indices.append(ii)
        if bc[2] != None:
            for ii in range(N, M*N - (2*N-1), N):
                boundaries.append(bc[2])
                indices.append(ii)
        if bc[3] != None:
            for ii in range(2*N-1, M*N - N, N):
                boundaries.append(bc[3])
                indices.append(ii)
        return indices, boundaries
    
if __name__ == '__main__':
        comm = MPI.COMM_WORLD
        app = Room(comm, 1/20.)
        solution = app.tempDistribution()
        data = comm.gather(solution,root=0)

        if comm.Get_rank() == 0:
            if 'saveData' in sys.argv:
                print 'Saves data to file...',
                with open('data.txt', 'w') as outfile:
                    # Converts data to JSON compatible format and dumps it in
                    # data.txt
                    json.dump([element.tolist() for element in data],
                              outfile, sort_keys = True, indent = 4,
                              ensure_ascii=False)
                print 'done!'

            if 'plot' in sys.argv:
                print 'Plots data...',
                Z = data[2]
                M = np.shape(Z)[0]
                N = np.shape(Z)[1]
                y = np.linspace(0,app.deltax * M,M)
                x = np.linspace(0,app.deltax * N,N)
                X, Y = np.meshgrid(x, y)

                plt.figure(1)
                minLim, maxLim, resolution = 4., 42., 20
                v = np.linspace(minLim, maxLim, resolution, endpoint=True)
                startpos=[[0.,0.],[1.,0],[2.,1]]
                endpos = [[1.,1.],[2.,2.],[3.,2.]]
                for ii in range(3):
                    Z = data[ii]
                    Z = np.flipud(Z)
                    M = np.shape(Z)[0]
                    N = np.shape(Z)[1]
                    y = np.linspace(startpos[ii][1], endpos[ii][1],M)
                    x = np.linspace(startpos[ii][0], endpos[ii][0],N)
                    X, Y = np.meshgrid(x, y)
                    plt.contourf(X, Y, Z, v)
                plt.colorbar(ticks=v)
                plt.xlabel('x')
                plt.ylabel('y')
                plt.show()
                print 'done!'
