#!/usr/bin/python3

""" B-Spline surface demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Importing general packages
# -------------------------------------------------------------------------------------------------------------------- #
import sys
import os
import time
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../')
from parablade.CAD_functions import BSplineSurface


# -------------------------------------------------------------------------------------------------------------------- #
# Define the B-Spline surface input parameters
# -------------------------------------------------------------------------------------------------------------------- #
# Array of control points
n_dim, n, m = 3, 5, 4
P = np.zeros((n_dim, n, m))

# First row
P[:, 0, 0] = [0.00, 1.00, 0.00]
P[:, 1, 0] = [1.00, 2.00, 0.00]
P[:, 2, 0] = [2.00, 2.50, 0.00]
P[:, 3, 0] = [3.00, 2.00, 0.00]
P[:, 4, 0] = [4.00, 1.00, 0.00]

# Second row
P[:, 0, 1] = [0.00, 1.00, 1.00]
P[:, 1, 1] = [1.00, 2.00, 1.00]
P[:, 2, 1] = [2.00, 2.50, 1.00]
P[:, 3, 1] = [3.00, 2.00, 1.00]
P[:, 4, 1] = [4.00, 1.00, 1.00]

# Third row
P[:, 0, 2] = [0.00, 1.00, 2.00]
P[:, 1, 2] = [1.00, 2.00, 2.00]
P[:, 2, 2] = [2.00, 2.50, 2.00]
P[:, 3, 2] = [3.00, 2.00, 2.00]
P[:, 4, 2] = [4.00, 1.00, 2.00]

# Fourth row
P[:, 0, 3] = [0.50, 1.00, 3.00]
P[:, 1, 3] = [1.00, 1.50, 3.00]
P[:, 2, 3] = [2.00, 2.00, 3.00]
P[:, 3, 3] = [3.00, 1.50, 3.00]
P[:, 4, 3] = [3.50, 1.00, 3.00]

# Maximum index of the control points (counting from zero)
n = np.shape(P)[1] - 1
m = np.shape(P)[2] - 1

# Define the order of the basis polynomials
# Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
# Set p = n (number of control points minus one) to obtain a Bezier
p = n
q = m

# Definition of the knot vectors (clamped spline)
# p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
# q+1 zeros, m minus q equispaced points between 0 and 1, and q+1 ones. In total s+1 points where s=m+q+1
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
V = np.concatenate((np.zeros(q), np.linspace(0, 1, m - q + 2), np.ones(q)))

# (u,v) parametrization. 1D arrays of (u,v) query points
Nu, Nv = 500, 500
u = np.linspace(0.00, 1.00, Nu)
v = np.linspace(0.00, 1.00, Nv)
[u,v] = np.meshgrid(u, v, indexing='xy')
u = u.flatten()
v = v.flatten()


# -------------------------------------------------------------------------------------------------------------------- #
# Create the B-Spline surface
# -------------------------------------------------------------------------------------------------------------------- #
t = time.time()
my_BSpline_surface = BSplineSurface(P, p, q, U, V)
S = my_BSpline_surface.get_BSplineSurface_value(u, v)
print('The elapsed time is %(my_time).3f seconds' % {'my_time': time.time() - t})


# -------------------------------------------------------------------------------------------------------------------- #
# Plot the B-Spline surface
# -------------------------------------------------------------------------------------------------------------------- #
options = {'point_cloud': 'no',
           'control_points': 'yes',
           'surface': 'yes',
           'surface_Nu': Nu,
           'surface_Nv': Nv}

my_BSpline_surface.plot_BSplineSurface(options)

plt.show()

