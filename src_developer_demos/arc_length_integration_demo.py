#!/usr/bin/python3

""" Arc length integration demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import numpy as np
import scipy.integrate as integrate
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../')
from parablade.CAD_functions import get_arc_length


# -------------------------------------------------------------------------------------------------------------------- #
# Function to be integrated
# -------------------------------------------------------------------------------------------------------------------- #
def my_helix(t):

    if np.isscalar(t):
        t = np.asarray([t])
    elif t.ndim == 1:
        t = t[np.newaxis]
    else:
        raise Exception('The parameter t must be an scalar or a one-dimensional array')

    x = 2/np.pi*np.cos(2*np.pi*t)
    y = 2/np.pi*np.sin(2*np.pi*t)
    z = 3*t
    C = np.concatenate((x,y,z), axis=0)
    return C


L_exact = 5.0
L_numeric = get_arc_length(my_helix, 0, 1)
print('The analytic arclength is', L_exact)
print('The numeric arclenth is', L_numeric)
print('The numerical error is', L_numeric-L_exact)



# -------------------------------------------------------------------------------------------------------------------- #
# Plot the parametric helix
# -------------------------------------------------------------------------------------------------------------------- #
# Prepare the plot
fig = plt.figure(figsize=(6,6))
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
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
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

# Draw line plot
t = np.linspace(0, 1, 500)
C = my_helix(t)
line, = ax.plot(C[0,:], C[1,:], C[2,:])
line.set_linewidth(1.25)
line.set_linestyle("-")
line.set_color("b")
line.set_marker(" ")
line.set_markersize(3.5)
line.set_markeredgewidth(1)
line.set_markeredgecolor("k")
line.set_markerfacecolor("w")
line.set_label('Circular helix')

# Create legend
ax.legend(ncol=1, loc='upper right', fontsize=13, edgecolor='k', framealpha=1.0)

# Set axes aspect ratio
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()
z_min, z_max = ax.get_zlim()
x_mid = (x_min + x_max) / 2
y_mid = (y_min + y_max) / 2
z_mid = (z_min + z_max) / 2
L = np.max((x_max - x_min, y_max - y_min, z_max - z_min)) / 2
ax.set_xlim3d(x_mid - 1.0 * L, x_mid + 1.0 * L)
ax.set_ylim3d(y_mid - 1.0 * L, y_mid + 1.0 * L)
ax.set_zlim3d(z_mid - 1.0 * L, z_mid + 1.0 * L)

# Adjust PAD
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Show the figure
plt.show()
