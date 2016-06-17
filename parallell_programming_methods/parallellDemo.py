 # -*- coding: utf-8 -*-
import os
import json
import numpy as np
import sys
"""
Usage: 
    1) Execute this script as is (``python parallellDemo.py''). The program
    will the pause while the plot is displayed. When closing the plot, the
    the data will load and be printed in the script.
    2) Use mpiexec (``mpiexec -np 3 python parallellDemo.py OPTIONS). The 
    processes will the print progress in the terminal. Note that OPTIONS can
    be defined as ``saveData'' to save the output data to a JSON file, or
    ``plot'', to visualize the results from the subprocesses themself. A
    combination of both can be used if the data should be both saved and
    plotted.
"""
print 'Solving problem...',
os.system('mpiexec -np 3 python parallellSolver.py saveData plot')
print 'complete!'

with open('data.txt', 'r') as infile:
    data = json.load(infile)
for ii in range(3):
    print '----------------- solution %s -----------------' % ii
    for jj in range(len(data[ii])):
        # Plotting data by row is necessary due to the size of the solution
        print data[ii][jj]
os.remove('data.txt')
