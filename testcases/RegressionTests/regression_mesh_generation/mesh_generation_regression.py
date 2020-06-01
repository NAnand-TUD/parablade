#!/usr/bin/python3

""" Regression rest to check if the mesh generation files are generated correctly """


#---------------------------------------------------------------------------------------------#
# Importing general packages
#---------------------------------------------------------------------------------------------#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
import sys
import os
import time


#---------------------------------------------------------------------------------------------#
# Generate the blade profile
#---------------------------------------------------------------------------------------------#
# To generate the mesh files run MakeBlade.py from terminal
# >> MakeBlade.py regression_blade.cfg


#---------------------------------------------------------------------------------------------#
# Read the mesh generation files
#---------------------------------------------------------------------------------------------#
# Hub surface
filename = os.getcwd() + '/output/mesh_files/' + 'hub.crv'
hub = np.loadtxt(filename, delimiter='\t').transpose()

# Shroud surface
filename = os.getcwd() + '/output/mesh_files/' + 'shroud.crv'
shroud = np.loadtxt(filename, delimiter='\t').transpose()

# Blade sections
filename = os.getcwd() + '/output/mesh_files/' + 'blade.crv'
blade_sections = np.loadtxt(filename, delimiter='\t', comments='#').transpose()


#---------------------------------------------------------------------------------------------#
# Print the mesh generation files
#---------------------------------------------------------------------------------------------#
# Create the figure
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
ax.set_xlabel('Axial direction', fontsize=11, color='k', labelpad=12)
ax.set_ylabel('Radial direction', fontsize=11, color='k', labelpad=12)

# Draw the hub surface
line, = ax.plot(hub[2, :], hub[0, :])
line.set_linewidth(1.25)
line.set_linestyle("-")
line.set_color("b")
line.set_marker(" ")
line.set_markersize(3.5)
line.set_markeredgewidth(1)
line.set_markeredgecolor("k")
line.set_markerfacecolor("w")
line.set_label(r'Hub surface')

# Draw the shroud surface
line, = ax.plot(shroud[2, :], shroud[0, :])
line.set_linewidth(1.25)
line.set_linestyle("-")
line.set_color("r")
line.set_marker(" ")
line.set_markersize(3.5)
line.set_markeredgewidth(1)
line.set_markeredgecolor("k")
line.set_markerfacecolor("w")
line.set_label(r'Shroud surface')

# Draw the blade sections
radial_coordinates = (blade_sections[0, :] ** 2 + blade_sections[1, :] ** 2) ** (1 / 2)
axial_coordinates = blade_sections[2, :]
line, = ax.plot(axial_coordinates, radial_coordinates)
line.set_linewidth(1.25)
line.set_linestyle(" ")
line.set_color("r")
line.set_marker("+")
line.set_markersize(3.5)
line.set_markeredgewidth(1)
line.set_markeredgecolor("k")
line.set_markerfacecolor("w")
line.set_label(r'Blade sections')

# Create legend
ax.legend(ncol=1, loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10, edgecolor='k', framealpha=1.0)

# Set the aspect ratio of the data
ax.set_aspect(1.0)

# Adjust PAD
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Show the figure
plt.show()


