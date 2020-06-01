#!/usr/bin/env python3

""" Convergence plot for the normal computation

    Compute the surface normals using different differentiation methods:

        1) Forward finite differences
        2) Central finite differences
        3) Complex step

    Plot the error of the differentiation as a function of the time step
    (assuming that the complex step with machine epsilon stepzsie is exact)

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
from common import printProgress
from config import ReadUserInput


# -------------------------------------------------------------------------------------------------------------------- #
# Preliminary definitions
# -------------------------------------------------------------------------------------------------------------------- #
# Create output directory if it does not exist
os.system("rm -rf figures")
os.mkdir(os.getcwd() + '/figures')

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

# Define a vector of stepsizes for the derivative computation
N_step = 100
step_sizes = np.logspace(-3, -14, N_step)


# -------------------------------------------------------------------------------------------------------------------- #
# Main computations
# -------------------------------------------------------------------------------------------------------------------- #
# Assume that the normal vectors computed using a machine epsilon complex step are exact
normals_exact = np.real(my_blade.get_surface_normals(u, v, method='complex_step', step=np.finfo(np.float64).eps))

# Initialize variables to print progress
total = len(step_sizes)
count = 0

# Compute the derivative error for each stepsize
error_CS = np.zeros((N_step,))
error_FFD = np.zeros((N_step,))
error_CFD = np.zeros((N_step,))
Normals_CS = np.zeros((3, N_step))
Normals_FFD = np.zeros((3, N_step))
Normals_CFD = np.zeros((3, N_step))
for h in step_sizes:

    # Compute the normal vectors with different methods using a range of step sizes
    normals_CS = np.real(my_blade.get_surface_normals(u, v, method='complex_step', step=h))
    normals_FFD = np.real(my_blade.get_surface_normals(u, v, method='forward_finite_differences', step=h))
    normals_CFD = np.real(my_blade.get_surface_normals(u, v, method='central_finite_differences', step=h))

    # Compute the two-norm of the error with respect to the "exact" solution
    error_CS[count] = np.sum((normals_CS - normals_exact)**2)**(1/2)
    error_FFD[count] = np.sum((normals_FFD - normals_exact) ** 2)**(1/2)
    error_CFD[count] = np.sum((normals_CFD - normals_exact) ** 2)**(1/2)

    # Append the sum of the (x,y,z) gradients over the (u,v) parameters for each stepsize
    # Summing the values for each (u,v) pair is an indicator to see if everything is correct for the entire blade
    Normals_CS[:, count] = np.sum(normals_CS** 2, axis=1)**(1/2)
    Normals_FFD[:, count] = np.sum(normals_FFD ** 2, axis=1)**(1/2)
    Normals_CFD[:, count] = np.sum(normals_CFD ** 2, axis=1)**(1/2)

    # Print progress bar
    printProgress(count, total)

    # Progress update
    count += 1


# -------------------------------------------------------------------------------------------------------------------- #
#  Plot the normal vector components
# -------------------------------------------------------------------------------------------------------------------- #
# Prepare figure and axes
fig, ax_list = plt.subplots(3, 1, sharex='col',figsize=(7.2, 6.2))
fontsize = 11
ax_list[0].set_title(r'2-norm of normal vector components', fontsize=fontsize + 1, pad=15)
x_label = r'$h$ - Stepsize'
y_labels = ['x component', 'y component', 'z component']

for i in range(3):
    # Select the axes on which to plot the results
    ax = ax_list[i]
    ax.set_xlim([10 ** -0, 10 ** -16])
    a = np.max(np.abs(np.concatenate((Normals_CS[i,:], Normals_FFD[i,:], Normals_CFD[i,:]))))-Normals_CS[i,-1]
    ax.set_ylim([Normals_CS[i,-1]+1.25*a, Normals_CS[i,-1]-1.25*a])
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_ylabel(y_labels[i], labelpad=10)
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0e'))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))
    if i == 2: ax.set_xlabel(x_label, fontsize=fontsize, color='k', labelpad=12)

    # Plot the error with respect to the "exact" result
    ax.plot(step_sizes, Normals_CS[i,:], linewidth=0.75, linestyle='-', color='k', label='CS')
    ax.plot(step_sizes, Normals_FFD[i,:], linewidth=0.75, linestyle='--', color='b', label='FFD')
    ax.plot(step_sizes, Normals_CFD[i,:], linewidth=0.75, linestyle='--', color='r', label='CFD')


# Prepare legend
ax_list[2].legend(edgecolor='k', fontsize=fontsize-1, loc='best', shadow = False, ncol = 3)

# Save the figure
fig.savefig('figures/normal_vector_components.png', bbox_inches='tight', dpi=200)


# ------------------------------------------------------------------------------------------------------------------ #
# Convergence plot
# ------------------------------------------------------------------------------------------------------------------ #
fig = plt.figure(figsize=(6, 5))
ax = fig.gca()
fontsize = 11
plt.title(r'Convergence plot for normal vector', pad=10, fontsize=fontsize+1)
ax.set_xlabel(r'$h$ - Stepsize', fontsize=fontsize, color='k', labelpad=12)
ax.set_ylabel(r'Two-norm of the error', fontsize=fontsize, color='k', labelpad=12)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.set_aspect('auto')
plt.yscale('log')
plt.xscale('log')
ax.set_xlim([10**-0, 10**-16])
# ax.set_ylim([10**-16, 10**1])

# Plot the error with respect to the "exact" result
line1, = ax.plot(step_sizes, error_CS, linewidth=0.75, linestyle='-', color='k', label='CS')
line2, = ax.plot(step_sizes, error_FFD, linewidth=0.75, linestyle='-', color='b', label='FFD')
line3, = ax.plot(step_sizes, error_CFD, linewidth=0.75, linestyle='-', color='r', label='CFD')

# Plot markers for the optimal stepsize (theoretical result, see Nocedal optimization textbook for details)
h_opt_FFD_double = (np.finfo(np.float64).eps)**(1/2)        # Forward finite difference optimal h
h_opt_CFD_double = (np.finfo(np.float64).eps)**(1/3)        # Central finite difference optimal h

normals_FFD = np.real(my_blade.get_surface_normals(u, v, method='forward_finite_differences', step=h_opt_FFD_double))
normals_CFD = np.real(my_blade.get_surface_normals(u, v, method='central_finite_differences', step=h_opt_CFD_double))
error_FFD_opt = np.sqrt(np.sum((normals_FFD - normals_exact) ** 2))
error_CFD_opt = np.sqrt(np.sum((normals_CFD - normals_exact) ** 2))

# Plot the points of theoretical optimum h for the finite difference methods
marker1, = ax.plot(h_opt_FFD_double,  error_FFD_opt,
                   linestyle=' ', color='b', marker='*', markerfacecolor='w', markersize=7, markeredgewidth=0.75,
                   label=r'FFD $h=\epsilon^{1/2}$')

marker2, = ax.plot(h_opt_CFD_double, error_CFD_opt,
                   linestyle=' ', color='r', marker='*', markerfacecolor='w', markersize=7, markeredgewidth=0.75,
                   label='CFD $h=\epsilon^{1/3}$')

# Prepare legend
ax.legend(edgecolor='k',loc='best',fontsize=fontsize-1)

# Adjust PAD
fig.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Save the figure
fig.savefig('figures/convergence_plot_normals.png', bbox_inches='tight', dpi=200)

plt.show()