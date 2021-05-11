#!/usr/bin/python3

""" Demonstration of the 2D blade parametrization based on connecting arcs """

#----------------------------------------------------------------------------------------------------------------------#
# Importing general packages
#----------------------------------------------------------------------------------------------------------------------#
import sys
import os
import time
import pdb
import numpy as np
import matplotlib.pyplot as plt


#----------------------------------------------------------------------------------------------------------------------#
# Importing user-defined packages
#----------------------------------------------------------------------------------------------------------------------#

from parablade.blade_2D_connecting_arcs import Blade2DConnectingArcs


#----------------------------------------------------------------------------------------------------------------------#
# Define the section design variables
#----------------------------------------------------------------------------------------------------------------------#
section_variables = {}
section_variables['stagger'] = -25
section_variables['theta_in'] = 37.7
section_variables['theta_out'] = -63.2
section_variables['wedge_in'] = 25 
section_variables['wedge_out'] = 5
section_variables['radius_in'] = 0.01
section_variables['radius_out'] = 0.010
section_variables['dist_1'] = 0.5
section_variables['dist_2'] = 0.30
section_variables['dist_3'] = 0.20
section_variables['dist_4'] = 0.20


# Convert standard-python scalars into singleton numpy arrays
for i in section_variables:
    section_variables[i] = np.asarray(section_variables[i])


#----------------------------------------------------------------------------------------------------------------------#
# Create the 2D section
#----------------------------------------------------------------------------------------------------------------------#
u = np.linspace(0.00, 1.00, 100)
my_blade = Blade2DConnectingArcs(section_variables)
my_blade.get_section_coordinates(u)


#----------------------------------------------------------------------------------------------------------------------#
# Plot the 2D section
#----------------------------------------------------------------------------------------------------------------------#
# Plot a single blade
my_blade.plot_blade(blade_section='yes', control_points='yes',
                    leading_edge_radius='yes', trailing_edge_radius='yes')

# Plot a cascade of blades
# my_blade.plot_blade_cascade()

# Plot the curvature distribution
# my_blade.plot_curvature_distribution()

# Show the figure
plt.show()


