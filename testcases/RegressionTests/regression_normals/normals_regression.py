#!/usr/bin/env python3

""" Regression rest for the surface normal vector computation

    Compute the surface normals using different differentiation methods:

        1) Forward finite differences
        2) Central finite differences
        3) Complex step

    Author: Roberto Agromayor
    Date: 28/06/2019
    Updated: 12/10/2019

"""


# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create output directory if it does not exist
if os.path.exists(os.getcwd() + '/figures/') is False:
    os.mkdir(os.getcwd() + '/figures/')


#----------------------------------------------------------------------------------------------------------------------#
# Import user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')
from blade_3D import Blade3D
from config import ReadUserInput


# -------------------------------------------------------------------------------------------------------------------- #
# Preliminary definitions
# -------------------------------------------------------------------------------------------------------------------- #
# Load the blade configuration file
IN_file = os.getcwd() + '/' + 'regression_blade.cfg'
IN = ReadUserInput(IN_file)

# (u,v) parametrization of the blade
h = 1e-3
Nu, Nv = 500, 50
u = np.linspace(0.00+h, 1.00-h, Nu)
v = np.linspace(0.00+h, 1.00-h, Nv)
[u, v] = np.meshgrid(u,v)
u = u.flatten()
v = v.flatten()
UV = np.concatenate((u[:, np.newaxis], v[:, np.newaxis]), axis=1)

# Create the blade
my_blade = Blade3D(IN, UV)
my_blade.make_blade()

# Get machine epsilon for double-precision floating-point arithmetics
eps = np.finfo(np.float64).eps


# -------------------------------------------------------------------------------------------------------------------- #
# Main computations
# -------------------------------------------------------------------------------------------------------------------- #
# Normal vector computation
normals_CS = my_blade.get_surface_normals(u, v, method='complex_step', step=1e-12)
normals_FFD = my_blade.get_surface_normals(u, v, method='forward_finite_differences', step=eps**(1/2))
normals_CFD = my_blade.get_surface_normals(u, v, method='central_finite_differences', step=eps**(1/3)/10)

# Assume that the normal vectors computed using a machine epsilon complex step are exact
normals_exact = my_blade.get_surface_normals(u, v, method='complex_step', step=eps)

# Two-norm of error
error_CS = np.real(np.sqrt(np.sum((normals_CS - normals_exact) ** 2))/Nu*Nv)
error_FFD = np.real(np.sqrt(np.sum((normals_FFD - normals_exact) ** 2))/Nu*Nv)
error_CFD = np.real(np.sqrt(np.sum((normals_CFD - normals_exact) ** 2))/Nu*Nv)

# Print the results
print('CS error:', error_CS, '\t FFD error:', error_FFD, '\t CFD error:', error_CFD)
