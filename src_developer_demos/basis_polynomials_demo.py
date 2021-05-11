#!/usr/bin/python3

""" B-Spline basis polynomials demonstration script """

# -------------------------------------------------------------------------------------------------------------------- #
# Importing general packages
# -------------------------------------------------------------------------------------------------------------------- #
import sys
import os
import time
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------------------------------------------- #
# Importing user-defined packages
# -------------------------------------------------------------------------------------------------------------------- #
sys.path.append(os.getcwd() + '/../')
from parablade.CAD_functions import BSplineCurve


# -------------------------------------------------------------------------------------------------------------------- #
# Compute the basis polynomials
# -------------------------------------------------------------------------------------------------------------------- #
# u-parametrization
Nu = 1000
u = np.linspace(0.00, 1, Nu)

# Array of control points (inmaterial for the computations)
P = np.zeros((1,5))

# Maximum index of the control points (counting from zero)
nn = np.shape(P)[1]
n = nn-1

# Define the order of the basis polynomials
# Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
# Set p = n (number of control points minus one) to obtain a Bezier
p = 3  # Set at most p = 3 (cubic B-Spline)

# Definition of the knot vectors (clamped spline)
# p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))


# Get the basis polynomials
my_BSpline = BSplineCurve(P, p, U)
C = my_BSpline.get_BSplineCurve_value(u)
N = np.real(my_BSpline.get_basis_polynomials(n, p, U, u))

# -------------------------------------------------------------------------------------------------------------------- #
# Plot the basis polynomials
# -------------------------------------------------------------------------------------------------------------------- #
# Create the figure
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111)
fontsize = 12
ax.set_xlabel('$u$ parameter', fontsize=fontsize, color='k', labelpad=12)
ax.set_ylabel('Basis polynomial value', fontsize=fontsize, color='k', labelpad=12)
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
# ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
# ax.set_xticks([])
# ax.set_yticks([])
# ax.axis('off')

for i in range(n+1):
    line, = ax.plot(u, N[i,:])
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    # line.set_color("k")
    line.set_marker(" ")
    line.set_markersize(3.5)
    line.set_markeredgewidth(1)
    line.set_markeredgecolor("k")
    line.set_markerfacecolor("w")
    line.set_label('index ' + str(i))

# Set the aspect ratio of the data
# ax.set_aspect(1.0)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Create legend
ax.legend(ncol=1, loc='best', fontsize=10, edgecolor='k', framealpha=1.0)

# Adjust pad
plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

# Show the figure
plt.show()


