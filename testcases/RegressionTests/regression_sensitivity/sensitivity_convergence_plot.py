#!/usr/bin/env python3

""" Convergence plot for the sensitivity computation

    Compute the sensitivity of the blade surface coordinates with respect to the design variables for different time
    steps using different differentiation methods:

        1) Forward finite differences
        2) Central finite differences
        3) Complex step

    Plot the error of the differentiation as a function of the time step
    (assuming that complex step with machine epsilon stepsize is exact)

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
N_step = 10
step_sizes = np.logspace(-2, -12, N_step)

# Chose what design variables to compute
my_names = my_blade.DVs_names
# my_names = my_blade.DVs_names_2D
# my_names = ['stagger']


# -------------------------------------------------------------------------------------------------------------------- #
# Main computations
# -------------------------------------------------------------------------------------------------------------------- #
# Compute the surface sensitivity for each design variable and stepsize
for key in my_names:
    my_numbers = range(len(my_blade.DVs_control_points[key]))
    for number in my_numbers:

        # Initialize variables to print progress
        print("Design variable:\t " + key + '_' + str(number))
        total = len(step_sizes)
        count = 0

        # Assume that the surface sensitivities computed using a machine epsilon complex step are exact
        grad_exact = my_blade.get_surface_sensitivity(u, v, method='complex_step', step=eps, variable=[key, number], display_progress='no')
        grad_exact = np.real(grad_exact[key + '_' + str(number)])
        Grad_exact = np.real(np.sum(grad_exact**2, axis=1)**(1/2))

        # Initialize arrays to store results
        Grad_CS = np.zeros((3, N_step))
        Grad_FFD = np.zeros((3, N_step))
        Grad_CFD = np.zeros((3, N_step))
        error_CS = np.zeros((N_step,))
        error_FFD = np.zeros((N_step,))
        error_CFD = np.zeros((N_step,))

        # Loop over the different step sizes
        for h in step_sizes:

            # Compute the surface sensitivity with different methods using a range of step sizes
            grad_CS = my_blade.get_surface_sensitivity(u, v, method='complex_step', step=h, variable=[key,number], display_progress='no')
            grad_FFD = my_blade.get_surface_sensitivity(u, v, method='forward_finite_differences', step=h, variable=[key,number], display_progress='no')
            grad_CFD = my_blade.get_surface_sensitivity(u, v, method='central_finite_differences', step=h, variable=[key,number], display_progress='no')

            # Convert to real numbers for the plotting
            grad_CS = np.real(grad_CS[key + '_' + str(number)])
            grad_FFD = np.real(grad_FFD[key + '_' + str(number)])
            grad_CFD = np.real(grad_CFD[key + '_' + str(number)])

            # Append the sum of the (x,y,z) gradients over the (u,v) parameters for each stepsize
            # Summing the values for each (u,v) pair is an indicator to see if everything is correct for the entire blade
            Grad_CS[:, count] = np.sum(grad_CS**2, axis=1)**(1/2)
            Grad_FFD[:, count] = np.sum(grad_FFD**2, axis=1)**(1/2)
            Grad_CFD[:, count] = np.sum(grad_CFD**2, axis=1)**(1/2)

            # Compute the error with respect to the "exact" derivative computation
            # The two-norm norm sums over the (x,y,z) and the (u,v) pairs (Frobenius Norm)
            # Summing the values for each (u,v) pair is an indicator to see if everything is correct for the entire blade
            error_CS[count] = np.sum((grad_CS - grad_exact)**2)**(1/2)
            error_FFD[count] = np.sum((grad_FFD - grad_exact)**2)**(1/2)
            error_CFD[count] = np.sum((grad_CFD - grad_exact)**2)**(1/2)

            # Print progress bar
            printProgress(count+1, total)

            # Progress update
            count += 1


        # ------------------------------------------------------------------------------------------------------------ #
        # Plot the gradient components
        # ------------------------------------------------------------------------------------------------------------ #
        # Prepare figure and axes
        fig, ax_list = plt.subplots(3, 1, sharex='col', figsize=(7.2, 6.2))
        fontsize = 11
        ax_list[0].set_title(r'Surface gradient with respect to: ' + key + '_' + str(number), fontsize=fontsize + 1, pad =15)
        x_label = r'$h$ - Stepsize'
        y_labels = ['$x$ component', '$y$ component', '$z$ component']

        for i in range(3):
            # Select the axes on which to plot the results
            ax = ax_list[i]
            ax.set_xlim([10 ** -0, 10 ** -16])
            a = np.max(np.abs(np.concatenate((Grad_CS[i,:], Grad_FFD[i,:], Grad_CFD[i,:]))))-Grad_exact[i]
            ax.set_ylim([Grad_exact[i]+1.25*a, Grad_exact[i]-1.25*a])
            ax.set_xscale('log')
            ax.set_yscale('linear')
            ax.set_ylabel(y_labels[i], labelpad=10)
            ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0e'))
            ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))
            if i == 2:
                ax.set_xlabel(x_label, fontsize=fontsize, color='k', labelpad=12)

            # Plot the error with respect to the "exact" result
            ax.plot(step_sizes, Grad_CS[i,:], linewidth=0.75, linestyle='-', color='k', label='CS')
            ax.plot(step_sizes, Grad_FFD[i,:], linewidth=0.75, linestyle='-.', color='b', label='FFD')
            ax.plot(step_sizes, Grad_CFD[i,:], linewidth=0.75, linestyle='--', color='r', label='CFD')

        # Prepare legend
        ax_list[2].legend(edgecolor='k', fontsize=fontsize-1, loc='lower left', shadow = False, ncol = 3)

        # Save the figure
        # fig.savefig('figures/gradient_components_'+key+'_'+str(number)+'.pdf', bbox_inches='tight')
        fig.savefig('figures/gradient_components_' + key + '_' + str(number) + '.png', bbox_inches='tight', dpi=200)


        # ------------------------------------------------------------------------------------------------------------ #
        # Convergence plot
        # ------------------------------------------------------------------------------------------------------------ #
        # Prepare figure and axes
        fig = plt.figure(figsize=(6, 5))
        ax = fig.gca()
        fontsize = 11
        plt.title(r'Convergence plot for: ' + key + '_' + str(number), pad=10, fontsize=fontsize + 1)
        ax.set_xlabel(r'$h$ - Stepsize', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel(r'Two-norm of the sensitivity error', fontsize=fontsize, color='k', labelpad=12)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.yscale('log')
        plt.xscale('log')
        ax.set_xlim([10 ** -0, 10 ** -16])
        ax.set_ylim([10**-17, 10**1])

        # Plot the error with respect to the "exact" result
        line1, = ax.plot(step_sizes, error_CS + eps, linewidth=0.75, linestyle='-', color='k', label='CS')
        line2, = ax.plot(step_sizes, error_FFD + 1.1 * eps, linewidth=0.75, linestyle='-.', color='b', label='FFD')
        line3, = ax.plot(step_sizes, error_CFD + 1.2 * eps, linewidth=0.75, linestyle='--', color='r', label='CFD')

        # Plot markers for the optimal stepsize (theoretical result, see Nocedal optimization textbook for details)
        h_opt_FFD_double = eps ** (1 / 2)  # Forward finite difference optimal h
        h_opt_CFD_double = eps ** (1 / 3)  # Central finite difference optimal h
        grad_FFD_opt = my_blade.get_surface_sensitivity(u, v, method='forward_finite_differences', step=h_opt_FFD_double, variable=[key,number], display_progress='no')
        grad_CFD_opt = my_blade.get_surface_sensitivity(u, v, method='central_finite_differences', step=h_opt_CFD_double, variable=[key,number], display_progress='no')
        error_FFD_opt = np.real(np.sum((grad_FFD_opt[key + '_' + str(number)] - grad_exact) ** 2)**(1/2))
        error_CFD_opt = np.real(np.sum((grad_CFD_opt[key + '_' + str(number)] - grad_exact) ** 2)**(1/2))

        # Plot the points of theoretical optimum h for the finite difference methods
        ax.plot(h_opt_FFD_double,  error_FFD_opt+eps,
                           linestyle=' ', color='b', marker='*', markerfacecolor='w', markersize=7, markeredgewidth=0.75,
                           label=r'FFD $h=\epsilon^{1/2}$')

        ax.plot(h_opt_CFD_double, error_CFD_opt+eps,
                           linestyle=' ', color='r', marker='*', markerfacecolor='w', markersize=7, markeredgewidth=0.75,
                           label='CFD $h=\epsilon^{1/3}$')

        # Prepare legend
        lgd = plt.legend(edgecolor='k', loc='best', fontsize=fontsize - 1)

        # Adjust PAD
        fig.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        # Save the figure
        # fig.savefig('figures/convergence_plot_' + key + '_' + str(number) + '.pdf', bbox_inches='tight')
        fig.savefig('figures/convergence_plot_' + key + '_' + str(number) + '.png', bbox_inches='tight', dpi=200)


# plt.show()
