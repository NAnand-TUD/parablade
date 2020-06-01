#!/usr/bin/python3

# """ Transfinite interpolation demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../')
from CAD_functions import BSplineCurve
from interpolation_functions import TransfiniteInterpolation


# -------------------------------------------------------------------------------------------------------------------- #
# West boundary parametrization
# -------------------------------------------------------------------------------------------------------------------- #
P = np.zeros((3, 4), dtype=complex)
P[:, 0] = [0.00, 0.00, 0.00]
P[:, 1] = [0.10, 0.33, +0.00]
P[:, 2] = [0.15, 0.66, -0.20]
P[:, 3] = [0.05, 1.00, 0.00]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c1 = BSplineCurve(P, p, U)
C1_func = c1.get_BSplineCurve_value
P1 = c1.P


# -------------------------------------------------------------------------------------------------------------------- #
# South boundary parametrization
# -------------------------------------------------------------------------------------------------------------------- #
P = np.zeros((3, 4), dtype=complex)
P[:, 0] = P1[:, 0]
P[:, 1] = [0.33, 0.00, -0.20]
P[:, 2] = [0.66, 0.10, 0.30]
P[:, 3] = [1.00, 0.20, 0.20]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c2 = BSplineCurve(P, p, U)
C2_func = c2.get_BSplineCurve_value
P2 = c2.P


# -------------------------------------------------------------------------------------------------------------------- #
# East boundary parametrization
# -------------------------------------------------------------------------------------------------------------------- #
P = np.zeros((3, 4), dtype=complex)
P[:, 0] = P2[:, -1]
P[:, 1] = [1.15, 0.33+0.25, 0.10]
P[:, 2] = [1.15, 0.66+0.25, 0.20]
P[:, 3] = [1.05, 1.00+0.25, 0.20]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c3 = BSplineCurve(P, p, U)
C3_func = c3.get_BSplineCurve_value
P3 = c3.P


# -------------------------------------------------------------------------------------------------------------------- #
# North boundary parametrization
# -------------------------------------------------------------------------------------------------------------------- #
P = np.zeros((3, 4), dtype=complex)
P[:, 0] = P1[:, -1]
P[:, 1] = [0.33, 1.15, 0.20]
P[:, 2] = [0.66, 1.15, 0.00]
P[:, 3] = P3[:, -1]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c4 = BSplineCurve(P, p, U)
C4_func = c4.get_BSplineCurve_value
P4 = c4.P


# -------------------------------------------------------------------------------------------------------------------- #
# Corners of the boundary
# -------------------------------------------------------------------------------------------------------------------- #
P12 = P1[:, 0]
P23 = P3[:, 0]
P34 = P3[:, -1]
P41 = P1[:, -1]


# -------------------------------------------------------------------------------------------------------------------- #
# Create the transfinite interpolator
# -------------------------------------------------------------------------------------------------------------------- #
transfinite_interpolator = TransfiniteInterpolation(C1_func, C2_func, C3_func, C4_func, P12, P23, P34, P41)


# -------------------------------------------------------------------------------------------------------------------- #
# Plot the domain boundary and control points
# -------------------------------------------------------------------------------------------------------------------- #
# Prepare the plot
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(azim=130, elev=15)
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
ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=10)
ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=10)
ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=10)
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

# Parameter to plot the boundary of the domain
u = np.linspace(0.00, 1.00, 1000)
v = np.linspace(0.00, 1.00, 1000)

# West boundary
line1, = ax.plot(np.real(C1_func(v)[0, :]), np.real(C1_func(v)[1, :]), np.real(C1_func(v)[2, :]))
line1.set_linewidth(2)
line1.set_linestyle("-")
line1.set_color("g")
line1.set_marker(" ")
line1.set_markersize(3.5)
line1.set_markeredgewidth(1)
line1.set_markeredgecolor("k")
line1.set_markerfacecolor("w")
line1.set_zorder(1)
line1.set_label(' ')

# South boundary
line2, = ax.plot(np.real(C2_func(u)[0, :]), np.real(C2_func(u)[1, :]), np.real(C2_func(v)[2, :]))
line2.set_linewidth(2)
line2.set_linestyle("-")
line2.set_color("r")
line2.set_marker(" ")
line2.set_markersize(3.5)
line2.set_markeredgewidth(1)
line2.set_markeredgecolor("k")
line2.set_markerfacecolor("w")
line2.set_zorder(1)
line2.set_label(' ')

# East boundary
line3, = ax.plot(np.real(C3_func(v)[0, :]), np.real(C3_func(v)[1, :]), np.real(C3_func(v)[2, :]))
line3.set_linewidth(2)
line3.set_linestyle("-")
line3.set_color("m")
line3.set_marker(" ")
line3.set_markersize(3.5)
line3.set_markeredgewidth(1)
line3.set_markeredgecolor("k")
line3.set_markerfacecolor("w")
line3.set_zorder(1)
line3.set_label(' ')

# North boundary
line4, = ax.plot(np.real(C4_func(u)[0, :]), np.real(C4_func(u)[1, :]), np.real(C4_func(v)[2, :]))
line4.set_linewidth(2)
line4.set_linestyle("-")
line4.set_color("b")
line4.set_marker(" ")
line4.set_markersize(3.5)
line4.set_markeredgewidth(1)
line4.set_markeredgecolor("k")
line4.set_markerfacecolor("w")
line4.set_zorder(1)
line4.set_label(' ')

# Plot the Coons patch as a wireframe
u = np.linspace(0.00, 1.00, 50)
v = np.linspace(0.00, 1.00, 50)
[u,v] = np.meshgrid(u, v)
u = u.flatten()
v = v.flatten()
C = transfinite_interpolator(u, v)
C = np.real(np.reshape(C, (3, 50, 50)))
ax.plot_wireframe(C[0, ...], C[1, ...], C[2, ...], alpha=0.25, linewidth=0.50, color='k', zorder=2)

# Adjust pad
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Show plot
plt.show()