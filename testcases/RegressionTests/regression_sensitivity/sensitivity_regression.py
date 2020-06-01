#!/usr/bin/env python3

""" Regression rest for the surface sensitivity computation

    Compute the sensitivity of the blade surface coordinates with respect to the design variables for different time
    steps using different differentiation methods:

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
h = 1.01e-3
Nu, Nv = 100, 50
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

## Chose what design variables to compute
my_names = my_blade.DVs_names
# my_names = my_blade.DVs_names_2D
# my_names = my_blade.DVs_names_meridional

# -------------------------------------------------------------------------------------------------------------------- #
# Main computations
# -------------------------------------------------------------------------------------------------------------------- #
# Compute the surface sensitivity for each design variable
print('{:>25} \t {:>20} \t {:>20} \t {:>20}'.format('Design variable', 'CS error', 'FFD error', 'CFD error'))
for key in my_names:
    my_numbers = range(len(my_blade.DVs_control_points[key]))
    for number in my_numbers:

        # Assume that the surface sensitivities computed using a machine epsilon complex step are exact
        grad_exact = my_blade.get_surface_sensitivity(u, v, method='complex_step', variable=[key,number], display_progress='no', step=eps)
        grad_exact = np.real(grad_exact[key + '_' + str(number)])

        # Compute the surface sensitivity with different methods using a suitable stepsize
        grad_CS  = my_blade.get_surface_sensitivity(u, v, method='complex_step',               variable=[key,number], display_progress='no', step=1e-12)
        grad_FFD = my_blade.get_surface_sensitivity(u, v, method='forward_finite_differences', variable=[key,number], display_progress='no', step=eps**(1/2))
        grad_CFD = my_blade.get_surface_sensitivity(u, v, method='central_finite_differences', variable=[key,number], display_progress='no', step=eps**(1/3))

        # Convert to real numbers
        grad_CS = np.real(grad_CS[key + '_' + str(number)])
        grad_FFD = np.real(grad_FFD[key + '_' + str(number)])
        grad_CFD = np.real(grad_CFD[key + '_' + str(number)])

        # Compute the error with respect to the "exact" derivative computation
        error_CS = np.sum((grad_CS - grad_exact) ** 2)**(1/2)
        error_FFD = np.sum((grad_FFD - grad_exact) ** 2)**(1/2)
        error_CFD = np.sum((grad_CFD - grad_exact) ** 2)**(1/2)

        # Print results
        print('{:>25} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format(key + '_' + str(number), error_CS, error_FFD, error_CFD))

