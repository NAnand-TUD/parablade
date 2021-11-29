#!/usr/bin/python3

""" Bicubic interpolation demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../parablade')
from parablade.interpolation_functions import BicubicInterpolation


# -------------------------------------------------------------------------------------------------------------------- #
# Function to be interpolated
# -------------------------------------------------------------------------------------------------------------------- #
def my_func(x1, x2):
    # z = 1 + x1**2 + x2**2
    z = 1 + np.sin(x1**2 + x2**2)
    # z = 1 + np.exp(x1+x2)
    # z = np.log(1 + x1**2 + x2**2)
    return z


# -------------------------------------------------------------------------------------------------------------------- #
# Original 2D grid and function values
# -------------------------------------------------------------------------------------------------------------------- #
a, b = 0.00, 2.00
Nx, Ny = 201, 101
x = np.linspace(a, b, Nx)
y = np.linspace(a, b, Ny)
[X,Y] = np.meshgrid(x, y, indexing='ij')
F = my_func(X, Y)


# -------------------------------------------------------------------------------------------------------------------- #
# Define a new grid for interpolation
# -------------------------------------------------------------------------------------------------------------------- #
a, b = 0.13, 1.89
nx, ny = 51, 51
xq = np.linspace(a, b, nx)
yq = np.linspace(a, b, ny)


# -------------------------------------------------------------------------------------------------------------------- #
# SciPy interpolation on a regular grid
# -------------------------------------------------------------------------------------------------------------------- #
time_start = time.time()
interpolant = RectBivariateSpline(x, y, F, kx=3, ky=3, s=0)
Fq = interpolant(xq, yq)
[Xq, Yq] = np.meshgrid(xq, yq, indexing='ij')
error = np.real(np.asscalar(np.sum((Fq - my_func(Xq, Yq))**2)**(1/2)/(nx*ny)))
print('SciPy interpolation: Error %(my_error).6e / Elapsed time %(my_time).6f seconds' % {'my_error': error, 'my_time': time.time() - time_start})


# -------------------------------------------------------------------------------------------------------------------- #
# ParaBlade bicubic interpolation
# -------------------------------------------------------------------------------------------------------------------- #
# Reshape the query points into 1D arrays
time_start = time.time()
xq = Xq.flatten()
yq = Yq.flatten()
f_interpolator = BicubicInterpolation(x, y, F)
fq = np.real(f_interpolator(xq,yq))

# Evaluate the interpolation error and running time
error = np.real(np.asscalar(np.sum((fq - my_func(xq, yq))**2)**(1/2)/(nx*ny)))
print('ParaBlade interpolation: Error %(my_error).6e / Elapsed time %(my_time).6f seconds' % {'my_error': error, 'my_time': time.time() - time_start})


# -------------------------------------------------------------------------------------------------------------------- #
# Plot the original and interpolated values
# -------------------------------------------------------------------------------------------------------------------- #
# Prepare the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(azim=120, elev=20)
ax.grid(False)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('k')
ax.yaxis.pane.set_edgecolor('k')
ax.zaxis.pane.set_edgecolor('k')
ax.xaxis.pane._alpha = 0.9
ax.yaxis.pane._alpha = 0.9
ax.zaxis.pane._alpha = 0.9
fontsize = 11
ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=12)
ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=12)
ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=12)
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
ax.xaxis.set_rotate_label(False)
ax.yaxis.set_rotate_label(False)
ax.zaxis.set_rotate_label(False)
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_zticks([])
# ax.axis('off')

# Plot the the original function as a surface
surf = ax.plot_surface(X, Y, F, alpha=0.5, color='b', linewidth=0, shade='no')

# # Plot the Scipy interpolated values as a cloud of points
# ax.plot(Xq.flatten(), Yq.flatten(), Fq.flatten() ,color='b', linestyle=' ', marker='o' ,markersize=1, markerfacecolor='b')

# Plot my bicubic interpolated values as a cloud of points
ax.plot(xq, yq, fq ,color='r', linestyle=' ', marker='o' ,markersize=1, markerfacecolor='r')

# Adjust pad
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Show plot
plt.show()


