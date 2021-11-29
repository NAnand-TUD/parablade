#!/usr/bin/python3

""" B-Spline curve demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Importing general packages
# -------------------------------------------------------------------------------------------------------------------- #
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../parablade')
from parablade.CAD_functions import BSplineCurve


# -------------------------------------------------------------------------------------------------------------------- #
# Create B-Spline curve object
# -------------------------------------------------------------------------------------------------------------------- #
# u-parametrization
u = np.linspace(0, 1, 1000)

# # Array of control points for a 2D case
# P = np.zeros((2,5), dtype=complex)
# P[:, 0] = [0.20, 0.50]
# P[:, 1] = [0.40, 0.70]
# P[:, 2] = [0.80, 0.60]
# P[:, 3] = [0.60, 0.20]
# P[:, 4] = [0.40, 0.20]

# Array of control points for a 3D case
P = np.zeros((3,5), dtype=complex)
P[:, 0] = [0.20, 0.50, 0.00]
P[:, 1] = [0.40, 0.70, 0.25]
P[:, 2] = [0.80, 0.60, 0.50]
P[:, 3] = [0.80, 0.40, 0.75]
P[:, 4] = [0.40, 0.20, 1.00]

# Maximum index of the control points (counting from zero)
nn = np.shape(P)[1]
n = nn-1

# Define the order of the basis polynomials
# Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
# Set p = n (number of control points minus one) to obtain a Bezier
p = 2

# Definition of the knot vectors (clamped spline)
# p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

# Create, evaluate, and plot the spline curve
my_BSpline = BSplineCurve(P, p, U)
C = my_BSpline.get_BSplineCurve_value(u)
options = {'line': 'yes', 'control_points': 'yes', 'order':0}
my_BSpline.plot_BSplineCurve(options)


plt.show()

