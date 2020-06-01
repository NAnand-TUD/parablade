###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

################################# FILE NAME: MakeBlade.py #####################################
#=============================================================================================#
# author: Roberto, Nitish Anand                                                               |
#    :PhD Candidates,                                                                         |
#    :Power and Propulsion, Energy Technology,                                                |
#    :NTNU, TU Delft,                                                                         |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#                                                                                             |
#=============================================================================================#

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import copy
import numpy as np


#----------------------------------------------------------------------------------------------------------------------#
# "Cluster mode" imports
#----------------------------------------------------------------------------------------------------------------------#
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except:
    pass


#----------------------------------------------------------------------------------------------------------------------#
# Import user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
from CAD_functions import *


# -------------------------------------------------------------------------------------------------------------------- #
# Define the 2D blade class (parametrization based on connecting arcs)
# -------------------------------------------------------------------------------------------------------------------- #
class Blade2DConnectingArcs:

    """ Create a 2D blade section object

    The parametrization is based on 4 points connected by 4 arcs

        - Leading edge (points 1-2)
        - Lower surface (points 2-3)
        - Trailing edge (points 3-4)
        - Upper surface (points 4-1)

    The blade section is G1 continuous at the connecting points

    Parameters
    ----------
    section_variables : dictionary of ndarrays with shape (n,)
    Dictionary containing the values of the blade section design variables:

        - Stagger angle
        - Inlet metal angle
        - Outlet metal angle
        - Leading edge wedge semi-angle
        - Trailing edge wedge semi-angle
        - Leading edge radius
        - Trailing edge radius
        - Tangent proportion at point 1
        - Tangent proportion at point 2
        - Tangent proportion at point 3
        - Tangent proportion at point 4

    References
    ----------
    TODO add publication information here

    """

    def __init__(self, section_variables):

        # Convert each singleton ndarray into a standard-python scalar
        section_variables = copy.deepcopy(section_variables)
        for i in section_variables:
            section_variables[i] = section_variables[i].item()

        # Declare input variables as instance variables
        self.stagger = section_variables["stagger"]
        self.theta_in = section_variables["theta_in"]
        self.theta_out = section_variables["theta_out"]
        self.wedge_in = section_variables["wedge_in"]
        self.wedge_out = section_variables["wedge_out"]
        self.radius_in = section_variables["radius_in"]
        self.radius_out = section_variables["radius_out"]
        self.dist_1 = section_variables["dist_1"]
        self.dist_2 = section_variables["dist_2"]
        self.dist_3 = section_variables["dist_3"]
        self.dist_4 = section_variables["dist_4"]

        # Declare additional variables as instance variables
        self.u = None
        self.leading_edge = None
        self.lower_side = None
        self.trailing_edge = None
        self.upper_side = None
        self.section_coordinates = None

        # Set the leading edge at (x,y)=(0,0)
        self.y_in = 0
        self.x_in = 0

        # Set the meridional chord equal to unity (normalized blade)
        self.chord = 1.00/np.cos(self.stagger)
        self.spacing = 0.75*self.chord

        # Compute the location of point 1
        self.theta1 = self.theta_in + self.wedge_in
        self.x1 = self.x_in + self.radius_in*(1-np.sin(self.theta1))
        self.y1 = self.y_in + self.radius_in*np.cos(self.theta1)

        # Compute the location of point 2
        self.theta2 = self.theta_in - self.wedge_in
        self.x2 = self.x_in + self.radius_in*(1+np.sin(self.theta2))
        self.y2 = self.y_in - self.radius_in*np.cos(self.theta2)

        # Compute the location of point 3
        self.theta3 = self.theta_out + self.wedge_out
        self.x3 = self.x_in + self.chord*np.cos(self.stagger) - self.radius_out*(1-np.sin(self.theta3))
        self.y3 = self.y_in + self.chord*np.sin(self.stagger) - self.radius_out *np.cos(self.theta3)

        # Compute the location of point 4
        self.theta4 = self.theta_out - self.wedge_out
        self.x4 = self.x_in + self.chord*np.cos(self.stagger) - self.radius_out*(1+np.sin(self.theta4))
        self.y4 = self.y_in + self.chord*np.sin(self.stagger) + self.radius_out*np.cos(self.theta4)

        # Create the arcs connecting points 1-2-3-4
        self.make_leading_edge()
        self.make_lower_side()
        self.make_trailing_edge()
        self.make_upper_side()


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the coodinates of the blade section
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_section_coordinates(self, u):

        """ Compute the coordinates of the blade section for input parameter u

        The u-parameter is divided such that

            - The leading edge is parametrized by ´u´ in the in interval [0.00, 0.25]
            - The lower surface is parametrized by ´u´ in the in interval [0.25, 0.50]
            - The trailing edge is parametrized by ´u´ in the in interval [0.50, 0.75]
            - The upper surface is parametrized by ´u´ in the in interval [0.75, 1.00]

        Parameters
        ----------
        u : ndarray with shape (N,)
            Array containing the u-parameter used to evaluate the B-Spline curve

        Returns
        -------
        section_coordinates : ndarray with shape (ndim, N)
            Array containing blade section coordinates
            The first dimension of ´section_coordinates´ spans the ´(x,y)´ coordinates
            The second dimension of ´section_coordinates´ spans the ´u´ parametrization sample points

        """

        # Sorting trick to retrieve order after concatenation
        my_order1 = np.argsort(u)
        my_order2 = np.argsort(my_order1)

        # Define the parameter for each arc (split the vector u in 4 pieces)
        # This is reordering, and we cannot reorder during matching. We retrieve the original order later
        u12 = np.sort((u[(u >= 0.00) & (u < 0.25)]-0.00)/(0.25-0.00))
        u23 = np.sort((u[(u >= 0.25) & (u < 0.50)]-0.25)/(0.50-0.25))
        u34 = np.sort((u[(u >= 0.50) & (u < 0.75)]-0.50)/(0.75-0.50))
        u41 = np.sort((u[(u >= 0.75) & (u <= 1.00)]-0.75)/(1.00-0.75))

        # Compute the coordinates of each arc
        coordinates_12 = self.leading_edge.get_NurbsCurve_value(u12)
        coordinates_23 = self.lower_side.get_NurbsCurve_value(u23)
        coordinates_34 = self.trailing_edge.get_NurbsCurve_value(u34)
        coordinates_41 = self.upper_side.get_NurbsCurve_value(u41)

        # Concatenate the arcs to obtain the blade section
        section_coordinates = np.concatenate((coordinates_12, coordinates_23, coordinates_34, coordinates_41), axis=1)

        # Retrieve the original parametrization order
        self.section_coordinates = section_coordinates[:, my_order2]

        return self.section_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the leading edge arc
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_leading_edge(self):

        """ Create the leading edge (circular arc parametrized as second order rational Bezier curve) """

        # Compute the control points
        P0 = np.array([self.x1, self.y1])                       # First control point
        P2 = np.array([self.x2, self.y2])                       # Third control point
        n = np.asarray([-(P2[1] - P0[1]), (P2[0] - P0[0])])     # Normal direction
        n = n / np.sqrt(n[0]**2+n[1]**2)                        # Unitary normal vector
        a = self.wedge_in                                       # Circular arc semi-angle
        R = self.radius_in                                      # Circular arc radius
        P1 = (P0 + P2) / 2 - R * np.cos(a) / np.tan(a) * n      # Second control point (see notes)
        P = np.array([P0, P1, P2]).transpose()                  # Matrix of control points

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Weight of the control points
        W = np.ones((n + 1,), dtype=complex)                    # Unitary weight
        W[1] = np.sin(a)                                        # Special weight for the second control point

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 2

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Create the leading edge NURBS curve object
        self.leading_edge = NurbsCurve(P, W, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the lower surface arc
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_lower_side(self):

        """ Create the lower surface (B-Spline of degree 3) """

        # Compute the control points
        P0 = np.array([self.x2, self.y2])
        P1 = np.array([self.x2, self.y2]) + self.dist_2 * self.chord * np.array(
            [np.cos(self.theta2), np.sin(self.theta2)])
        P2 = np.array([self.x3, self.y3]) - self.dist_3 * self.chord * np.array(
            [np.cos(self.theta3), np.sin(self.theta3)])
        P3 = np.array([self.x3, self.y3])
        P = np.array([P0, P1, P2, P3]).transpose()              # Matrix of control points

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Weight of the control points
        W = np.ones((n + 1,), dtype=complex)                    # Unitary weight

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Create the lower surface NURBS curve object
        self.lower_side = NurbsCurve(P, W, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the trailing edge arc
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_trailing_edge(self):

        """ Create the trailing edge (circular arc parametrized as second order rational Bezier curve) """

        # Compute the control points
        P0 = np.array([self.x3, self.y3])                       # First control point
        P2 = np.array([self.x4, self.y4])                       # Third control point
        n = np.asarray([-(P2[1] - P0[1]), (P2[0] - P0[0])])     # Normal direction
        n = n / np.sqrt(n[0]**2+n[1]**2)                        # Unitary normal vector
        a = self.wedge_out                                      # Circular arc semi-angle
        R = self.radius_out                                     # Circular arc radius
        P1 = (P0 + P2) / 2 - R * np.cos(a) / np.tan(a) * n      # Second control point (see notes)
        P = np.array([P0, P1, P2]).transpose()                  # Matrix of control points

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Weight of the control points
        W = np.ones((n + 1,), dtype=complex)                    # Unitary weight
        W[1] = np.sin(a)                                        # Special weight for the second control point

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 2

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Create the trailing edge NURBS curve object
        self.trailing_edge = NurbsCurve(P, W, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the upper surface arc
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_upper_side(self):

        """ Create the upper surface (B-Spline of degree 3) """

        # Compute the control points
        P0 = np.array([self.x4, self.y4])
        P1 = np.array([self.x4, self.y4]) - self.dist_4 * self.chord * np.array([np.cos(self.theta4), np.sin(self.theta4)])
        P2 = np.array([self.x1, self.y1]) + self.dist_1 * self.chord * np.array([np.cos(self.theta1), np.sin(self.theta1)])
        P3 = np.array([self.x1, self.y1])
        P = np.array([P0, P1, P2, P3]).transpose()              # Matrix of control points

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Weight of the control points
        W = np.ones((n + 1,))                                   # Unitary weight

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Create the upper surface NURBS curve object
        self.upper_side = NurbsCurve(P, W, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute section curvature
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_section_curvature(self, u, h=1e-5):

        """ Compute the curvature of the blade section """

        # Compute curvature by central finite differences
        dr = (self.get_section_coordinates(u+h) - self.get_section_coordinates(u))/h
        ddr = (self.get_section_coordinates(u+h) - 2*self.get_section_coordinates(u) + self.get_section_coordinates(u-h))/(h**2)
        section_curvature = np.abs(np.cross(ddr, dr, axisa=0, axisb=0)) / (np.sum(dr ** 2, axis=0)) ** (3 / 2)

        return section_curvature


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plot the blade section
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_blade(self, fig=None, ax=None,
                   blade_section='yes', control_points='yes',
                   leading_edge_radius='no', trailing_edge_radius='no'):

        """ Plot of the 2D blade section (control appearance using the options dictionary) """

        # Create the figure
        if fig is None:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)

        fontsize = 12
        ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.axis('off')
        h = 0
        u_plot = np.linspace(0.00+h, 1-h, 300)

        # Plot the blade surface
        if blade_section == 'yes':

            # Plot leading edge
            line, = ax.plot(np.real(self.leading_edge.get_NurbsCurve_value(u_plot)[0, :]),
                            np.real(self.leading_edge.get_NurbsCurve_value(u_plot)[1, :]))
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot lower_surface
            line, = ax.plot(np.real(self.lower_side.get_NurbsCurve_value(u_plot)[0, :]),
                            np.real(self.lower_side.get_NurbsCurve_value(u_plot)[1, :]))
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot trailing edge
            line, = ax.plot(np.real(self.trailing_edge.get_NurbsCurve_value(u_plot)[0, :]),
                            np.real(self.trailing_edge.get_NurbsCurve_value(u_plot)[1, :]))
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot leading edge line
            line, = ax.plot(np.real(self.upper_side.get_NurbsCurve_value(u_plot)[0, :]),
                            np.real(self.upper_side.get_NurbsCurve_value(u_plot)[1, :]))
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')


        # Plot the blade section control polyline
        if control_points == 'yes':

            # Plot leading edge control points
            line, = ax.plot(np.real(self.leading_edge.P[0, :]), np.real(self.leading_edge.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot lower surface control points
            line, = ax.plot(np.real(self.lower_side.P[0, :]), np.real(self.lower_side.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot trailing edge control points
            line, = ax.plot(np.real(self.trailing_edge.P[0, :]), np.real(self.trailing_edge.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot upper surface control points
            line, = ax.plot(np.real(self.upper_side.P[0, :]), np.real(self.upper_side.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')


            # Plot connecting points
            line, = ax.plot([self.x1, self.x2, self.x3, self.x4], [self.y1, self.y2, self.y3, self.y4])
            line.set_linewidth(0.75)
            line.set_linestyle(" ")
            line.set_color("k")
            line.set_marker("o")
            line.set_markersize(4)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')


        if leading_edge_radius == 'yes':

            # Get the circle coordinates
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.radius_in
            y_c = self.y_in
            x_circle = x_c + self.radius_in * np.cos(theta)
            y_circle = y_c + self.radius_in * np.sin(theta)

            # Draw the circle
            line, = ax.plot(x_circle, y_circle)
            line.set_linewidth(0.75)
            line.set_linestyle("-")
            line.set_color("b")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')


        if trailing_edge_radius == 'yes':

            # Get the circle coordinates
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.chord * np.cos(self.stagger) - self.radius_out
            y_c = self.y_in + self.chord * np.sin(self.stagger)
            x_circle = x_c + self.radius_out * np.cos(theta)
            y_circle = y_c + self.radius_out * np.sin(theta)

            # Draw the circle
            line, = ax.plot(x_circle, y_circle)
            line.set_linewidth(0.75)
            line.set_linestyle("-")
            line.set_color("b")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Set the aspect ratio of the data
        ax.set_aspect(1.0)

        # # Set the aspect ratio of the figure
        # ratio = 1.00
        # x1, x2 = ax.get_xlim()
        # y1, y2 = ax.get_ylim()
        # ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

        # # Adjust pad
        # plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        return fig, ax


    def plot_blade_cascade(self, fig=None, ax=None):

        """ Plot the blade cascade """

        # Create the figure
        if fig is None:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)

        fontsize = 12
        ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.axis('off')
        u_plot = np.linspace(0.00, 1.00, 500)

        for k in range(3):

            # Draw blade surface
            line, = ax.plot(np.real(self.get_section_coordinates(u_plot)[0, :]),
                            np.real(self.get_section_coordinates(u_plot)[1, :])+ self.spacing * (k - 1))
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Set the aspect ratio of the data
        ax.set_aspect(1.0)

        # # Set the aspect ratio of the figure
        # ratio = 1.00
        # x1, x2 = ax.get_xlim()
        # y1, y2 = ax.get_ylim()
        # ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

        # # Adjust pad
        # plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        return fig, ax


    def plot_curvature_distribution(self, fig=None, ax=None):

        """ Plot the curvature distribution of the 2D parametrization """

        # Create the figure
        if fig is None:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)

        fontsize = 12
        ax.set_xlabel('$u$ parameter', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel('$\kappa$ - Curvature', fontsize=fontsize, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.axis('off')

        # Parameter vector for the plotting
        h = 1e-5
        hh = h + h**2
        u12 = np.linspace(0.00 + hh, 0.25 - hh, 200)
        u23 = np.linspace(0.25 + hh, 0.50 - hh, 200)
        u34 = np.linspace(0.50 + hh, 0.75 - hh, 200)
        u41 = np.linspace(0.75 + hh, 1.00 - hh, 200)

        # Leading edge curvature
        section_curvature = np.real(self.get_section_curvature(u12, h))
        line, = ax.plot(u12, section_curvature)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("b")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Leading edge')

        # Lower surface curvature
        section_curvature = np.real(self.get_section_curvature(u23, h))
        line, = ax.plot(u23, section_curvature)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("r")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Lower surface')

        # Trailing edge curvature
        section_curvature = np.real(self.get_section_curvature(u34, h))
        line, = ax.plot(u34, section_curvature)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("g")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Trailing edge')

        # Upper surface curvature
        section_curvature = np.real(self.get_section_curvature(u41, h))
        line, = ax.plot(u41, section_curvature)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Upper surface')

        # # Set the aspect ratio of the data
        # ax.set_aspect(1.0)

        # Set the aspect ratio of the figure
        ratio = 1.00
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        # Create legend
        ax.legend()

        return fig, ax

