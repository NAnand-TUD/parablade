#!/usr/bin/env python3

""" Test the effect of each design variable on a 2D blade section

    Run this script to visually check the influence of each design variable on a 2D blade
    This script can be used to gain insight on the influence of each variable on the geometry of the blade or to check
    that the parametrization is behaving as expected

    This script can be used for blades using both the CONNECTING_ARCS or the CAMBER_THICKNESS parametrizations

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
# Plot the influence of the design variables
# -------------------------------------------------------------------------------------------------------------------- #
# Create output directory if it does not exist
os.system("rm -rf figures")
os.mkdir(os.getcwd() + '/figures')

# Compute the baseline geometry
IN_file = os.getcwd() + '/' + 'LS89_2D_connecting_arcs.cfg'
# IN_file = os.getcwd() + '/' + 'LS89_2D_camber_thickness.cfg'
IN = ReadUserInput(IN_file)
my_blade = Blade3D(IN)
my_blade.make_blade()
IN = my_blade.IN
coordinates_baseline = copy.deepcopy(np.real(my_blade.surface_coordinates))

# Select what variables are analyzed
# variable_names = my_blade.DVs_names_2D
variable_names = my_blade.DVs_names
# variable_names = ['stagger']

# Compute the geometry when each design variable is increased and decreased
for k in variable_names:
    for i in range(len(IN[k])):

        # Increase the design variables by an arbitrary amount
        IN_bis = copy.deepcopy(IN)
        if k in ['x_leading', 'y_leading', 'z_leading', 'x_trailing', 'z_trailing',
                 'x_hub', 'z_hub', 'x_shroud', 'z_shroud']:
            IN_bis[k][i] = IN[k][i] + 0.005
        elif k in ['theta_in', 'theta_out', 'stagger']:
            IN_bis[k][i] = IN[k][i] + 5
        elif k in ['wedge_in', 'wedge_out']:
            IN_bis[k][i] = IN[k][i] + 5
        elif k in ['radius_in', 'radius_out']:
            IN_bis[k][i] = IN[k][i] * 1.5
        elif k in ['dist_in', 'dist_out', 'dist_1', 'dist_2', 'dist_3', 'dist_4']:
            IN_bis[k][i] = IN_bis[k][i] + 0.20
        elif k in ['thickness_upper_1', 'thickness_upper_2', 'thickness_upper_3',
                   'thickness_upper_4', 'thickness_upper_5', 'thickness_upper_6',
                   'thickness_lower_1', 'thickness_lower_2', 'thickness_lower_3',
                   'thickness_lower_4', 'thickness_lower_5', 'thickness_lower_6']:
            IN_bis[k][i] = IN[k][i] + 0.05

        # Compute the new geometry
        my_blade.update_DVs_control_points(IN_bis)
        my_blade.make_blade()
        coordinates_increment = copy.deepcopy(np.real(my_blade.surface_coordinates))

        # Decrease the design variables by an arbitrary amount
        IN_bis = copy.deepcopy(IN)
        if k in ['x_leading', 'y_leading', 'z_leading', 'x_trailing', 'z_trailing',
                 'x_hub', 'z_hub', 'x_shroud', 'z_shroud']:
            IN_bis[k][i] = IN[k][i] - 0.005
        elif k in ['theta_in', 'theta_out', 'stagger']:
            IN_bis[k][i] = IN[k][i] - 5
        elif k in ['wedge_in', 'wedge_out']:
            IN_bis[k][i] = IN[k][i] - 5
        elif k in ['radius_in', 'radius_out']:
            IN_bis[k][i] = IN[k][i] / 1.5
        elif k in ['dist_in', 'dist_out', 'dist_1', 'dist_2', 'dist_3', 'dist_4']:
            IN_bis[k][i] = IN_bis[k][i] - 0.20
        elif k in ['thickness_upper_1', 'thickness_upper_2', 'thickness_upper_3',
                   'thickness_upper_4', 'thickness_upper_5', 'thickness_upper_6',
                   'thickness_lower_1', 'thickness_lower_2', 'thickness_lower_3',
                   'thickness_lower_4', 'thickness_lower_5', 'thickness_lower_6']:
            IN_bis[k][i] = IN[k][i] - 0.05

        # Compute the new geometry
        my_blade.update_DVs_control_points(IN_bis)
        my_blade.make_blade()
        coordinates_decrement = copy.deepcopy(np.real(my_blade.surface_coordinates))

        # Plot the results for each design variable
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(111)
        fontsize = 12
        plt.xlabel('$x$ axis (mm)', fontsize=fontsize, color='k', labelpad=12)
        plt.ylabel('$y$ axis (mm)', fontsize=fontsize, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        ax.axis('equal')

        # Plot the blades
        ax.plot(1e3*coordinates_baseline[0, :], 1e3*coordinates_baseline[1, :], linestyle='-', linewidth=0.75,  color='k', markeredgewidth=0.75, marker=' ', markersize=3,
                markerfacecolor='w', label='Baseline geometry')
        ax.plot(1e3*coordinates_increment[0, :], 1e3*coordinates_increment[1, :], linestyle='-', color='r', linewidth=0.75,  markeredgewidth=0.75, marker=' ', markersize=3,
                markerfacecolor='w', label='Increase '+k)
        ax.plot(1e3*coordinates_decrement[0, :], 1e3*coordinates_decrement[1, :], linestyle='-', color='b', linewidth=0.75, markeredgewidth=0.75, marker=' ', markersize=3,
                markerfacecolor='w', label='Decrease '+k)

        # Make the legend
        plt.legend(loc='best', ncol=1, edgecolor='k', fontsize=10)

        # Adjust PAD
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        # Save the figures
        fig.savefig('figures/influence_'+k+'.png', bbox_inches='tight', dpi=200)

