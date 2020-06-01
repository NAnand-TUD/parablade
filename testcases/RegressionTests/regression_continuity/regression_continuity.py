#!/usr/bin/env python3

""" Regression rest to check if the the parametrization is continuous at the edges:

    Connecting arcs parametrization:

        Edge 1:     u = 0.00|1.00       Upper surface / Leading edge
        Edge 2:     u = 0.25            Leading edge / Lower surface
        Edge 3:     u = 0.50            Lower surface / Trailing edge
        Edge 4:     u = 0.75            Trailing edge / Upper surface

    Camber thickness parametrization

        Edge 1:     u = 0.00|1.00       Upper surface / Lower surface at leading edge
        Edge 2:     u = 0.50            Upper surface / Lower surface at trailing edge


    The script tests the sensitivity in terms of

    1) Surface coordinates
    2) Surface normals
    3) Surface sensitivity with respect to each design variable


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

# (u,v) parametrization of the blade section
h = 1e-3
u_lower_approximation = [1.00-h, 0.25-h, 0.50-h, 0.75-h, 1.00-h]
u_upper_approximation = [0.00+h, 0.25+h, 0.50+h, 0.75+h, 0.00+h]
N = 5

# Choose the variables to investigate
blade_lower = Blade3D(IN)
my_names = blade_lower.DVs_names_2D

# Initialize lists to store the results for each edge (each value of my_u)
continuity_coordinates = []
continuity_normals = []
continuity_sensitivity = []


# -------------------------------------------------------------------------------------------------------------------- #
# Main computations
# -------------------------------------------------------------------------------------------------------------------- #
# Compute the continuity error for each edge
for i in range(N):

    # Define the (u,v) parametrization
    Nu, Nv = N, 100
    u1 = np.asarray(u_lower_approximation[i])
    u2 = np.asarray(u_upper_approximation[i])
    v = np.linspace(0.001, 0.999, 50)
    [u1, v1] = np.meshgrid(u1, v)
    [u2, v2] = np.meshgrid(u2, v)
    u1, v1 = u1.flatten(), v1.flatten()
    u2, v2 = u2.flatten(), v2.flatten()
    UV1 = np.concatenate((u1[:, np.newaxis]-h, v1[:, np.newaxis]), axis=1)    # Lower approximation
    UV2 = np.concatenate((u2[:, np.newaxis]+h, v2[:, np.newaxis]), axis=1)    # Upper approximation

    # Make a blade with the left approximation (u,v) parametrization
    blade_lower = Blade3D(IN, UV1)
    blade_lower.make_blade()

    # Make a blade with the right approximation (u,v) parametrization
    blade_upper = Blade3D(IN, UV2)
    blade_upper.make_blade()

    # Compute the difference in coordinates between lower and upper approximations to the edge
    # 2-norm of the XYZ coordinates and v-parametrization
    continuity_coordinates.append(np.real(np.sum((blade_lower.surface_coordinates - blade_upper.surface_coordinates)**2) ** (1/2) / Nv))

    # Compute the difference in normal vectors between lower and upper approximations to the edge
    # 2-norm of the XYZ coordinates and v-parametrization
    continuity_normals.append(np.sum((blade_lower.surface_normals - blade_upper.surface_normals)**2) ** (1/2) / Nv)

    # Compute the difference in sensitivity between lower and upper approximations to the edge (each design variable)
    # 2-norm of the XYZ coordinates and v-parametrization
    sensitivity_dict = {}
    for key in my_names:
        my_numbers = range(len(blade_lower.DVs_control_points[key]))
        for number in my_numbers:
            grad_1 = blade_lower.surface_sensitivity[key + '_' + str(number)]
            grad_2 = blade_upper.surface_sensitivity[key + '_' + str(number)]
            sensitivity_dict[key + '_' + str(number)] = np.sum((grad_1 - grad_2)**2)**(1/2)/Nv

    continuity_sensitivity.append(sensitivity_dict)



# -------------------------------------------------------------------------------------------------------------------- #
# Print the results
# -------------------------------------------------------------------------------------------------------------------- #
# Print headers
print('{:>20} \t {:>20} \t {:>20} \t {:>20} \t {:>20} \t {:>20}'.format('Sensitivity type',
                                                                        'u = 0.00',
                                                                        'u = 0.25',
                                                                        'u = 0.50',
                                                                        'u = 0.75',
                                                                        'u = 1.00'))

# Print coordinate mismatch
print('{:>20} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format('Coordinates',
                                                                                       continuity_coordinates[0],
                                                                                       continuity_coordinates[1],
                                                                                       continuity_coordinates[2],
                                                                                       continuity_coordinates[3],
                                                                                       continuity_coordinates[4]))
# Print normal vectors mismatch
print('{:>20} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format('Normals',
                                                                                       continuity_normals[0],
                                                                                       continuity_normals[1],
                                                                                       continuity_normals[2],
                                                                                       continuity_normals[3],
                                                                                       continuity_normals[4]))

# Print sensitivity mismatch
for key in my_names:
    my_numbers = range(len(blade_lower.DVs_control_points[key]))
    for number in my_numbers:
        print('{:>20} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format(key + '_' + str(number),
                                              continuity_sensitivity[0][key + '_' + str(number)],
                                              continuity_sensitivity[1][key + '_' + str(number)],
                                              continuity_sensitivity[2][key + '_' + str(number)],
                                              continuity_sensitivity[3][key + '_' + str(number)],
                                              continuity_sensitivity[4][key + '_' + str(number)]))
