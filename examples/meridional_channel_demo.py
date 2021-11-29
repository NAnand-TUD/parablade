#!/usr/bin/python3

""" Meridional channel parametrization demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Importing general packages
# -------------------------------------------------------------------------------------------------------------------- #
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../parablade')
from parablade.CAD_functions import BSplineCurve
from parablade.interpolation_functions import TransfiniteInterpolation

# -------------------------------------------------------------------------------------------------------------------- #
# Meridional channel parametrization
# -------------------------------------------------------------------------------------------------------------------- #

# West boundary parametrization
P = np.zeros((2,4), dtype=complex)
P[:, 0] = [0.00, 0.00]
P[:, 1] = [0.10, 0.33]
P[:, 2] = [0.15, 0.66]
P[:, 3] = [0.05, 1.00]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c1 = BSplineCurve(P, p, U)
C1_func = c1.get_BSplineCurve_value
P1 = c1.P

# South boundary parametrization
P = np.zeros((2,4), dtype=complex)
P[:, 0] = P1[:, 0]
P[:, 1] = [0.33, 0.00]
P[:, 2] = [0.66, 0.10]
P[:, 3] = [1.00, 0.20]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c2 = BSplineCurve(P, p, U)
C2_func = c2.get_BSplineCurve_value
P2 = c2.P

# East boundary parametrization
P = np.zeros((2,4), dtype=complex)
P[:, 0] = P2[:, -1]
P[:, 1] = [1.15, 0.33+0.25]
P[:, 2] = [1.15, 0.66+0.25]
P[:, 3] = [1.05, 1.00+0.25]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c3 = BSplineCurve(P, p, U)
C3_func = c3.get_BSplineCurve_value
P3 = c3.P

# North boundary parametrization
P = np.zeros((2,4), dtype=complex)
P[:, 0] = P1[:, -1]
P[:, 1] = [0.33, 1.15]
P[:, 2] = [0.66, 1.15]
P[:, 3] = P3[:, -1]
nn = np.shape(P)[1]
n = nn-1
p = n
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
c4 = BSplineCurve(P, p, U)
C4_func = c4.get_BSplineCurve_value
P4 = c4.P

# Boundary corners
P12 = P1[:, 0]
P23 = P3[:, 0]
P34 = P3[:, -1]
P41 = P1[:, -1]

# Transfinite interpolator
transfinite_interpolator = TransfiniteInterpolation(C1_func, C2_func, C3_func, C4_func, P12, P23, P34, P41)


# -------------------------------------------------------------------------------------------------------------------- #
# Plot the domain boundary and control points
# -------------------------------------------------------------------------------------------------------------------- #
# Prepare the figure
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
fontsize = 12
ax.set_xlabel('x - axis', fontsize=fontsize, color='k', labelpad=12)
ax.set_ylabel('y - axis', fontsize=fontsize, color='k', labelpad=12)
ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)

# Plot the boundary of the domain
u = np.linspace(0.00, 1.00, 1000)
v = np.linspace(0.00, 1.00, 1000)
line1, = ax.plot(C1_func(v)[0, :], C1_func(v)[1, :])
line2, = ax.plot(C2_func(u)[0, :], C2_func(u)[1, :])
line3, = ax.plot(C3_func(v)[0, :], C3_func(v)[1, :])
line4, = ax.plot(C4_func(u)[0, :], C4_func(u)[1, :])
lines = (line1, line2, line3, line4)
for line in lines:
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    line.set_color("k")
    line.set_marker(" ")
    line.set_markersize(3.5)
    line.set_markeredgewidth(1)
    line.set_markeredgecolor("k")
    line.set_markerfacecolor("w")
    line.set_label(' ')

# Plot the boundary control polylines
# Modify the control point coordinates to make plot more pretty
p1 = P1 + np.asarray([[+0.00, +0.05, +0.02, +0.00], [+0.00, +0.00, +0.00, +0.00]])
p2 = P2 + np.asarray([[+0.00, +0.00, +0.00, +0.00], [+0.00, -0.03, +0.06, +0.00]])
p3 = P3 + np.asarray([[+0.00, +0.03, +0.03, +0.00], [+0.00, +0.00, +0.00, +0.00]])
p4 = P4 + np.asarray([[+0.00, +0.00, +0.00, +0.00], [+0.00, +0.01, -0.03, +0.00]])

line1, = ax.plot(p1[0, :], p1[1, :])
line2, = ax.plot(p2[0, :], p2[1, :])
line3, = ax.plot(p3[0, :], p3[1, :])
line4, = ax.plot(p4[0, :], p4[1, :])
lines = (line1, line2, line3, line4)
for line in lines:
    line.set_linewidth(0.50)
    line.set_linestyle("-")
    line.set_color("k")
    line.set_marker("o")
    line.set_markersize(3.5)
    line.set_markeredgewidth(0.75)
    line.set_markeredgecolor("k")
    line.set_markerfacecolor("w")
    line.set_label(' ')


# # Plot the u-direction lines
# u = np.linspace(0.00, 1.00, 1000)
# v = np.linspace(0.00, 1.00, 8)
# for v_value in v:
#     C = transfinite_interpolator(u, v_value)
#     line, = ax.plot(C[0, :], C[1, :])
#     line.set_linewidth(0.50)
#     line.set_linestyle("-")
#     line.set_color("k")
#     line.set_marker(" ")
#     line.set_markersize(3.5)
#     line.set_markeredgewidth(1)
#     line.set_markeredgecolor("r")
#     line.set_markerfacecolor("w")
#     line.set_label(' ')
#
# # Plot the v-direction lines
# u = np.linspace(0.00, 1.00, 8)
# v = np.linspace(0.00, 1.00, 1000)
# for u_value in u:
#     C = transfinite_interpolator(u_value, v)
#     line, = ax.plot(C[0, :], C[1, :])
#     line.set_linewidth(0.50)
#     line.set_linestyle("-")
#     line.set_color("k")
#     line.set_marker(" ")
#     line.set_markersize(3.5)
#     line.set_markeredgewidth(1)
#     line.set_markeredgecolor("r")
#     line.set_markerfacecolor("w")
#     line.set_label(' ')


# Set the aspect ratio of the data
ax.set_aspect(1.0)

# Adjust pad
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Hide axes
plt.axis('off')

# Save the figure
# plt.savefig('figures/meridional_channel.svg', bbox_inches='tight')
# plt.savefig('figures/meridional_channel.pdf', bbox_inches='tight')

# Show the figure
plt.show()


