# Parallellized solver for the 2D-laplace equation in a cartesian coordinate system
This is the third project in the course in 'Advanced Numerical Algorithms'
(FMNN25) and concerns the implementation of parallellized solver for a 2D-heat
problem on an irregular grid consisting of three rooms [[1]]. The approach makes use
the Python Message Passing Interface (MPI) and shows how parallellized solvers
can be implemented and run.

The script `parallellSolver.py' can be from the terminal launched using the command,

  ``mpiexec -np 3 python parallellSolver.py OPTIONS'',

where OPTIONS are command line arguments set to 'plot' to visualize the solution,
'saveData' to save the data to a file as a JSON, a combination of both arguments
or nothing at all.

The script `parallellDemo.py' can be run from the terminal as a regular python script
and will then solve the problem and print the data.

The software makes use of the python MPI package [[2]] to parallellize the solver. It
is straight forward to install on a Linux or OS X operating system [[3]], but slightly harder
on windows. A good guide can be found here [[4]].

The solver is not yet general enough to encompass all types of 2D-heat problems
capable of being represented by many rectangular sub-grids. In order to achieve a
more general, version the following will be required:

* Initiate the rooms from some interface (the demo script or a GUI) and not in the
MPI processes themselves.
* Pass an argument with boundary normals (currently rank/room-specific).
* Handle corners and check for incompatibilities in sub-grids.

[1]: http://www.maths.lth.se/na/courses/FMNN25/media/material/projectDDPython.pdf
[2]: https://mpi4py.scipy.org/docs/usrman/tutorial.html
[3]: https://pypi.python.org/pypi/mpi4py
[4]: http://www.maths.lth.se/na/courses/FMNN25/media/material/MPI.pdf
