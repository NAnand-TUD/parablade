#!/usr/bin/python3

""" Demonstration of the 2D blade parametrization based on connecting arcs """

#----------------------------------------------------------------------------------------------------------------------#
# Importing general packages
#----------------------------------------------------------------------------------------------------------------------#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


#----------------------------------------------------------------------------------------------------------------------#
# Importing user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
sys.path.append(os.getcwd() + '/../parablade')
from parablade.blade_2D_connecting_arcs import Blade2DConnectingArcs


#----------------------------------------------------------------------------------------------------------------------#
# Define the section design variables
#----------------------------------------------------------------------------------------------------------------------#
section_variables = {}
section_variables['stagger'] = -45 * np.pi / 180
section_variables['theta_in'] = +20 * np.pi / 180
section_variables['theta_out'] = -70 * np.pi / 180
section_variables['wedge_in'] = 25 * np.pi / 180
section_variables['wedge_out'] = 6 * np.pi / 180
section_variables['radius_in'] = 0.100
section_variables['radius_out'] = 0.020
section_variables['dist_1'] = 0.35
section_variables['dist_2'] = 0.30
section_variables['dist_3'] = 0.40
section_variables['dist_4'] = 0.70

# Convert standard-python scalars into singleton numpy arrays
for i in section_variables:
    section_variables[i] = np.asarray(section_variables[i])


#----------------------------------------------------------------------------------------------------------------------#
# Create the 2D section
#----------------------------------------------------------------------------------------------------------------------#
u = np.linspace(0.00, 1.00, 1000)
my_blade = Blade2DConnectingArcs(section_variables)
my_blade.get_section_coordinates(u)


#----------------------------------------------------------------------------------------------------------------------#
# Plot the 2D section
#----------------------------------------------------------------------------------------------------------------------#
# Plot a single blade
my_blade.plot_blade(blade_section='yes', control_points='yes',
                    leading_edge_radius='no', trailing_edge_radius='no')

# Plot a cascade of blades
my_blade.plot_blade_cascade()

# Plot the curvature distribution
my_blade.plot_curvature_distribution()

# Show the figure
plt.show()


