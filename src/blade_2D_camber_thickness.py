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
# Define the 2D blade class (parametrization based on camberline and thickness)
# -------------------------------------------------------------------------------------------------------------------- #
class Blade2DCamberThickness:

    """ Create a 2D blade section object

    The parametrization is based on two arcs (upper surface and lower surface)
    The arcs are constructed imposing a thickness distribution normal to the camberline
    The thickness distributions of the upper and lower surfaces are independent
    The connection between upper and lower surfaces at the leading/trailing edges is G2 continuous

    Parameters
    ----------
    section_variables : dictionary of ndarrays with shape (n,)
    Dictionary containing the values of the blade section design variables:

        - Stagger angle
        - Inlet metal angle
        - Outlet metal angle
        - Leading edge radius
        - Trailing edge radius
        - Camberline leading edge tangent proportion
        - Camberline trailing edge tangent proportion
        - Upper surface thickness distribution control points
        - Lower surface thickness distribution control points

    References
    ----------
    TODO add publication information here

    """

    def __init__(self, section_variables):

        # Convert each singleton ndarray into a standard-python scalar
        section_variables = copy.deepcopy(section_variables)
        for i in section_variables:
            section_variables[i] = section_variables[i].item()

        # Load the blade section variables
        self.stagger = section_variables["stagger"]
        self.theta_in = section_variables["theta_in"]
        self.theta_out = section_variables["theta_out"]
        self.radius_in = section_variables["radius_in"]
        self.radius_out = section_variables["radius_out"]
        self.dist_in = section_variables["dist_in"]
        self.dist_out = section_variables["dist_out"]

        self.thickness_upper = [self.radius_in,
                                section_variables['thickness_upper_1'],
                                section_variables['thickness_upper_2'],
                                section_variables['thickness_upper_3'],
                                section_variables['thickness_upper_4'],
                                section_variables['thickness_upper_5'],
                                section_variables['thickness_upper_6'],
                                self.radius_out]

        self.thickness_lower = [self.radius_in,
                                section_variables['thickness_lower_1'],
                                section_variables['thickness_lower_2'],
                                section_variables['thickness_lower_3'],
                                section_variables['thickness_lower_4'],
                                section_variables['thickness_lower_5'],
                                section_variables['thickness_lower_6'],
                                self.radius_out]

        # Declare additional variables as instance variables
        self.u = None
        self.u_sample_points = None
        self.camberline_sample_points = None
        self.camberline = None
        self.upper_thickness_distribution = None
        self.lower_thickness_distribution = None
        self.upper_side = None
        self.lower_side_BSpline = None
        self.section_coordinates = None

        # Set the leading edge at (x,y)=(0,0)
        self.y_in = 0
        self.x_in = 0

        # Set the meridional chord equal to unity (normalized blade)
        self.chord = 1.00 / np.cos(self.stagger)
        self.spacing = 0.75 * self.chord

        # Create camberline
        self.make_camberline()

        # Create thickness distribution
        self.make_upper_thickness_distribution()
        self.make_lower_thickness_distribution()

        # Get the sampling points to impose the thickness distribution
        self.N_sample = 10
        self.u_sample_points = self.get_sampling_points(sampling_mode='uniform')

        # Create blade surfaces
        self.make_upper_side()
        self.make_lower_side()


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the coodinates of the blade section
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_section_coordinates(self, u):

        """ Compute the coordinates of the blade section for input parameter u

        The u-parameter is divided such that

            - The upper surface is parametrized by ´u´ in the in interval [0.00, 0.50]
            - The lower surface is parametrized by ´u´ in the in interval [0.50, 1.00]

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
        u_lower = np.sort((u[(u >= 0.00) & (u < 0.50)] - 0.00) / (0.50 - 0.00))
        u_upper = np.sort((u[(u >= 0.50) & (u <= 1.00)] - 0.50) / (1.00 - 0.50))

        # Compute the coordinates of each arc
        upper_surface_coordinates = self.get_upper_side_coordinates(u_upper)
        lower_surface_coordinates = self.get_lower_side_coordinates(u_lower)

        # Concatenate the arcs to obtain the blade section
        section_coordinates = np.concatenate((lower_surface_coordinates, upper_surface_coordinates), axis=1)

        # Retrieve the original parametrization order
        self.section_coordinates = section_coordinates[:, my_order2]

        return self.section_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the camberline
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_camberline(self):

        """ Create the camberline curve (B-Spline of degree 3) """

        # Control point coordintes
        x0 = self.x_in
        y0 = self.y_in
        x1 = self.x_in + self.dist_in * np.cos(self.theta_in)
        y1 = self.y_in + self.dist_in * np.sin(self.theta_in)
        x3 = self.x_in + self.chord * np.cos(self.stagger)
        y3 = self.x_in + self.chord * np.sin(self.stagger)
        x2 = x3 - self.dist_out * np.cos(self.theta_out)
        y2 = y3 - self.dist_out * np.sin(self.theta_out)

        # Array of control points
        P = np.asarray([[x0, y0], [x1, y1], [x2, y2], [x3, y3]]).transpose()

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Get the B-Spline curve object
        self.camberline = BSplineCurve(P, p, U)


    def get_camberline_coordinates(self, u):

        """ Evaluate and return the coordinates of the camberline for input u """
        camberline_coordinates = self.camberline.get_BSplineCurve_value(u)

        return camberline_coordinates


    def get_camberline_slope(self, u):

        """ Evaluate and return the slope of the camberline for input u """
        camberline_slope = self.camberline.get_BSplineCurve_derivative(u, order=1)

        return camberline_slope


    def get_camberline_normal(self, u):

        """ Evaluate and return the unitary vector normal to the camber line """
        dr_sample = self.camberline.get_BSplineCurve_derivative(u, order=1)
        tangent = dr_sample / np.sum(dr_sample ** 2, axis=0) ** (1 / 2)
        camberline_normal = np.asarray([-tangent[1, :], tangent[0, :]])

        return camberline_normal


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the upper and lower thickness distributions
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_upper_thickness_distribution(self):

        """ Create the upper thickness distribution (B-Spline of degree 3) """

        # Array of control points
        P = np.asarray([self.thickness_upper])

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Get the B-Spline curve object
        self.upper_thickness_distribution = BSplineCurve(P, p, U)


    def make_lower_thickness_distribution(self):

        """ Create the thickness distribution (B-Spline of degree 3) """

        # Array of control points
        P = np.asarray([self.thickness_lower])

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Get the B-Spline curve object
        self.lower_thickness_distribution = BSplineCurve(P, p, U)


    def get_upper_thickness_distribution_values(self, u):

        """ Evaluate and return the coordinates of the camberline for input u """
        upper_thickness_distribution = self.upper_thickness_distribution.get_BSplineCurve_value(u)

        return upper_thickness_distribution


    def get_lower_thickness_distribution_values(self, u):

        """ Evaluate and return the coordinates of the camberline for input u """
        lower_thickness_distribution = self.lower_thickness_distribution.get_BSplineCurve_value(u)

        return lower_thickness_distribution


    def get_sampling_points(self, sampling_mode='uniform'):

        """ Distribute the sampling points along the camberline

        Parameters
        ----------
        sampling_mode : string
            Type of clustering used to sample the thickness distribution

        Returns
        -------
        u_sample : ndarray with shape (N_sample,)
            Array containing the values of the u-parameter used to sample the thickness distributions

        Notes
        -----
        The following options are available for the variable ´sampling_mode´:
            - Uniform spacing of the sampling points
            - Cosine-sigmoid distribution that clusters the sample points towards the leading and trailing edges
            - Sigmoid distribution that clusters the sample points towards the leading and trailing edges

        """

        if sampling_mode == 'uniform':
            u_sample = np.linspace(0, 1, self.N_sample)

        elif sampling_mode == 'cosine':
            u = np.linspace(0, 1, self.N_sample)
            u_sample = np.cos(np.pi / 2 * (1 - u)) ** 2

        elif sampling_mode == 'cluster':
            u = np.linspace(1e-6, 1 - 1e-6, self.N_sample)
            beta = 3 / 2  # Set beta > 0  | Larger values of beta increase clustering
            u_sample = 1 / (1 + (u / (1 - u)) ** (-beta))

        else:
            raise Exception("Choose a valid option for the sampling distribution: 'uniform', 'cosine', 'cluster'")

        return u_sample


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the upper surface
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_upper_side(self):

        # Number of sampling points used to impose the thickness distribution
        N = self.N_sample
        n = N + 2 - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Get the camberline coordinates, thickness, and normal vector at the sampling points
        camberline_sample = self.get_camberline_coordinates(self.u_sample_points)
        self.camberline_sample_points = camberline_sample
        normal_sample = self.get_camberline_normal(self.u_sample_points)  # Positive sign for upper surface
        thickness_sample = self.get_upper_thickness_distribution_values(self.u_sample_points)

        # Get the upper surface set of control points
        P = camberline_sample + normal_sample * thickness_sample

        # Collapse the first and last control points into the camberline
        P[:, 0] = camberline_sample[:, 0]
        P[:, -1] = camberline_sample[:, -1]

        # Get additional control points to ensure metal angles an G2 continuity
        P1 = self.get_start_G2_control_point(P, p, U, normal_sample[:, 0], 1 / self.radius_in)
        P2 = self.get_end_G2_control_point(P, p, U, normal_sample[:, -1], 1 / self.radius_out)

        # Concatenate all the control points
        P = np.concatenate((P[:, 0, np.newaxis], P1, P[:, 1:-1], P2, P[:, -1, np.newaxis]), axis=1)

        # Reverse the order of the control points such that the surface is parametrized counter-clockwise for u in [0,1]
        P = P[:, ::-1]

        # Get the B-Spline curve object
        self.upper_side = BSplineCurve(P, p, U)


    def get_upper_side_coordinates(self, u):

        """ Evaluate and return the coordinates of the upper_surface for input u """
        upper_side_coordinates = self.upper_side.get_BSplineCurve_value(u)

        return upper_side_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the lower surface
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_lower_side(self):

        # Number of sampling points used to impose the thickness distribution
        N = self.N_sample
        n = N + 2 - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 3

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Get the camberline coordinates, thickness, and normal vector at the sampling points
        camberline_sample = self.get_camberline_coordinates(self.u_sample_points)
        normal_sample = self.get_camberline_normal(self.u_sample_points) * (-1)  # Negative sign for lower surface
        thickness_sample = self.get_lower_thickness_distribution_values(self.u_sample_points)

        # Get the upper surface array initial set of control points
        P = camberline_sample + normal_sample * thickness_sample

        # Collapse the first and last control points into the camberline
        P[:, 0] = camberline_sample[:, 0]
        P[:, -1] = camberline_sample[:, -1]

        # Get additional control points to ensure metal angles an G2 continuity
        P1 = self.get_start_G2_control_point(P, p, U, normal_sample[:, 0], 1 / self.radius_in)
        P2 = self.get_end_G2_control_point(P, p, U, normal_sample[:, -1], 1 / self.radius_out)

        # Concatenate all the control points
        P = np.concatenate((P[:, 0, np.newaxis], P1, P[:, 1:-1], P2, P[:, -1, np.newaxis]), axis=1)

        # Get the B-Spline curve object
        self.lower_side_BSpline = BSplineCurve(P, p, U)


    def get_lower_side_coordinates(self, u):

        """ Evaluate and return the coordinates of the lower_surface for input u """
        lower_side_coordinates = self.lower_side_BSpline.get_BSplineCurve_value(u)

        return lower_side_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Define functions to ensure G2 continuity
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_start_G2_control_point(self, P, p, U, normal, curvature):

        """ Get the location of the additional control point that imposes the desired curvature at ´u=0´

        Compute the coordinates of the control point ´P1´ such that the B-Spline defined
        by the set of control points ´P´ has the specified curvature at ´u=0´

        Parameters
        ----------
        P : ndarray with shape (ndim=2, n)
            Array containing the coordinates of the original set of control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

        p : int
            Degree of the B-Spline curve

        U : ndarray with shape (r+1=n+p+2,)
            Array with the knot vector of the modified B-Spline curve

        normal: ndarray with shape (ndim=2,)
            Vector normal to the B-Spline curve at the start-point ´u=0´

        curvature: scalar
            Desired curvature at the the start-point ´u=0´
            The curvature is the inverse of the radius of curvature: \kappa = 1/radius

        Returns
        -------
        P1: ndarray with shape (ndim=2,)
            Array containing the coordinates of the extra control point that ensures G2 continuity

        """

        # First derivative
        # a_0 = -p / U[p + 1]
        a_1 = +p / U[p + 1]

        # Second derivative
        k = p * (p - 1) / U[p + 1]
        # b_0 = 1 / U[p + 1] * k
        # b_1 = -(1 / U[p + 1] + 1 / U[p + 2]) * k
        b_2 = 1 / U[p + 2] * k

        # Analytic expression to impose the radius of curvature of a B-Spline curve
        P0 = P[:, 0]
        P2 = P[:, 1]
        # dist_02 = np.linalg.norm(P2 - P0)
        dist_02 = np.sum((P2-P0)**2)**(1/2)
        alpha = np.arccos(np.dot(P2 - P0, normal) / dist_02)
        dist_01 = np.sqrt(b_2 / a_1 ** 2 / curvature * dist_02 * np.sin(alpha))
        P1 = P0 + normal * dist_01
        P1 = P1[:, np.newaxis]

        return P1


    def get_end_G2_control_point(self, P, p, U, normal, curvature):

        """ Get the location of the additional control point that imposes the desired curvature at ´u=1´

        Compute the coordinates of the control point ´P1´ such that the B-Spline defined
        by the set of control points ´P´ has the specified curvature at ´u=1´

        Parameters
        ----------
        P : ndarray with shape (ndim=2, n)
            Array containing the coordinates of the original set of control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

        p : int
            Degree of the B-Spline curve

        U : ndarray with shape (r+1=n+p+2,)
            Array with the knot vector of the modified B-Spline curve

        normal: ndarray with shape (ndim=2,)
            Vector normal to the B-Spline curve at the end-point ´u=1´

        curvature: scalar
            Desired curvature at the the end-point ´u=1´
            The curvature is the inverse of the radius of curvature: \kappa = 1/radius

        Returns
        -------
        P1: ndarray with shape (ndim=2,)
            Array containing the coordinates of the extra control point that ensures G2 continuity

        """

        # First derivative
        m = len(U) - 1
        n = m - p - 1
        # a_0 = +p / (1 - U[n])
        a_1 = -p / (1 - U[n])

        # Second derivative
        k = p * (p - 1) / (1 - U[n])
        # b_0 = 1 / (1 - U[n]) * k
        # b_1 = -(1 / (1 - U[n]) + 1 / (1 - U[n - 1])) * k
        b_2 = 1 / (1 - U[n - 1]) * k

        # Analytic expression to impose the radius of curvature of a B-Spline curve
        P0 = P[:, -1]
        P2 = P[:, -2]
        # dist_02 = np.linalg.norm(P2 - P0)
        dist_02 = np.sum((P2 - P0) ** 2) ** (1 / 2)
        alpha = np.arccos(np.dot(P2 - P0, normal) / dist_02)
        dist_01 = np.sqrt(b_2 / a_1 ** 2 / curvature * dist_02 * np.sin(alpha))
        P1 = P0 + normal * dist_01
        P1 = P1[:, np.newaxis]

        return P1


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute section curvature
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_section_curvature(self, u, h=1e-5):

        """ Compute the curvature of the blade section by central finite differences

        Parameters
        ----------
        u : ndarray with shape (Nu,)
            Array containing the u-parameter used to evaluate the section curvature

        h : scalar
            Step size used for the central differences
            Too large values lead to large truncation errors
            Too small values lead to large rounding errors

        Returns
        -------
        section_curvature : ndarray with shape (Nu, )
            Array containing the curvature distribution of the blade section

        """

        # Compute curvature by central finite differences
        dr = (self.get_section_coordinates(u+h) - self.get_section_coordinates(u))/h
        ddr = (self.get_section_coordinates(u+h) - 2*self.get_section_coordinates(u) + self.get_section_coordinates(u-h))/(h**2)
        section_curvature = np.abs(np.cross(ddr, dr, axisa=0, axisb=0)) / (np.sum(dr ** 2, axis=0)) ** (3 / 2)

        return section_curvature


    def check_analytic_curvature(self):

        """ Compute the curvature at the leading/trailing edges analytically and compare it with the input value """

        # Leading and trailing edge u-values
        u = np.linspace(0, 1, 2)

        # Upper surface curvature
        dr = self.upper_side.get_BSplineCurve_derivative(u, order=1)
        ddr = self.upper_side.get_BSplineCurve_derivative(u, order=2)
        upper_curvature = np.real(np.abs(np.cross(ddr, dr, axisa=0, axisb=0)) / (np.sum(dr ** 2, axis=0)) ** (3 / 2))

        # Lower surface curvature
        dr = self.lower_side_BSpline.get_BSplineCurve_derivative(u, order=1)
        ddr = self.lower_side_BSpline.get_BSplineCurve_derivative(u, order=2)
        lower_curvature = np.real(np.abs(np.cross(ddr, dr, axisa=0, axisb=0)) / (np.sum(dr ** 2, axis=0)) ** (3 / 2))

        # Print curvature mismatch
        print('\tCheck the the curvature given as a design variable matches with the output curvature:')
        print('{:>20} \t {:>20}    \t {:>20}    \t {:>20}'.format('Point', 'Design variable', 'Upper surface analytic', 'Lower surface analytic'))
        print('{:>20} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format('Start', 1/self.radius_in,  upper_curvature[1], lower_curvature[0]))
        print('{:>20} \t {:>20.5e} \t {:>20.5e} \t {:>20.5e}'.format('End',   1/self.radius_out, upper_curvature[0], lower_curvature[1]))


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_blade_section(self, fig=None, ax=None,
                           upper_side='yes', upper_side_control_points='no',
                           lower_side='yes', lower_side_control_points='no',
                           camberline='no', camberline_control_points='no', camberline_sample_points='no',
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

        # Plot camber line
        if camberline == 'yes':
            camberline_coordinates = np.real(self.get_camberline_coordinates(u_plot))
            line, = ax.plot(camberline_coordinates[0, :], camberline_coordinates[1, :])
            line.set_linewidth(0.75)
            line.set_linestyle("--")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw upper side
        if upper_side == 'yes':
            upper_surface_coordinates = np.real(self.get_upper_side_coordinates(u_plot))
            line, = ax.plot(upper_surface_coordinates[0, :], upper_surface_coordinates[1, :])
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw lower side
        if lower_side == 'yes':
            lower_surface_coordinates = np.real(self.get_lower_side_coordinates(u_plot))
            line, = ax.plot(lower_surface_coordinates[0, :], lower_surface_coordinates[1, :])
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw upper side control points
        if upper_side_control_points == 'yes':
            line, = ax.plot(np.real(self.upper_side.P[0, :]), np.real(self.upper_side.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label('o')

        # Draw lower side control points
        if lower_side_control_points == 'yes':
            line, = ax.plot(np.real(self.lower_side_BSpline.P[0, :]), np.real(self.lower_side_BSpline.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw the sample points on the camber line
        if camberline_sample_points == 'yes':
            line, = ax.plot(np.real(self.camberline_sample_points[0, :]), np.real(self.camberline_sample_points[1, :]))
            line.set_linewidth(1.0)
            line.set_linestyle(" ")
            line.set_color("k")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw the control points defining the camber line
        if camberline_control_points == 'yes':
            line, = ax.plot(np.real(self.camberline.P[1, :]), np.real(self.camberline.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw the osculating center at the leading edge
        if leading_edge_radius == 'yes':
            # Get the circle coordinates
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.radius_in * np.cos(self.theta_in)
            y_c = self.y_in + self.radius_in * np.sin(self.theta_in)
            x_circle = x_c + self.radius_in * np.cos(theta)
            y_circle = y_c + self.radius_in * np.sin(theta)

            # Draw the circle
            line, = ax.plot(x_circle, y_circle)
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("b")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label(' ')

        # Draw the osculating circle at the trailing edge
        if trailing_edge_radius == 'yes':
            # Get the circle coordinates
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.chord * np.cos(self.stagger) - self.radius_out * np.cos(self.theta_out)
            y_c = self.y_in + self.chord * np.sin(self.stagger) - self.radius_out * np.sin(self.theta_out)
            x_circle = x_c + self.radius_out * np.cos(theta)
            y_circle = y_c + self.radius_out * np.sin(theta)

            # Draw the circle
            line, = ax.plot(x_circle, y_circle)
            line.set_linewidth(1.25)
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

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        return fig, ax

    def plot_blade_cascade(self, fig=None, ax=None):

        """ Plot the cascade of blades """

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

            # Draw upper surface
            section_coordinates = np.real(self.get_section_coordinates(u_plot))
            line, = ax.plot(section_coordinates[0, :], section_coordinates[1, :] + self.spacing * (k - 1))
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

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        # Hide axes
        plt.axis('off')

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
        u_upper = np.linspace(0.00 + hh, 0.50 - hh, 1000)
        u_lower = np.linspace(0.50 + hh, 1.00 - hh, 1000)

        # Plot upper surface curvature
        line, = ax.plot(2*u_upper, np.real(self.get_section_curvature(u_upper)))
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("b")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Upper surface')

        # Plot lower surface curvature
        line, = ax.plot(2*u_upper[::-1], np.real(self.get_section_curvature(u_lower)))
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("r")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Lower surface')

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


    def plot_thickness_distribution(self, fig=None, ax=None):

        """ Plot of the 2D thickness distribution """

        # Create the figure
        if fig is None:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)

        fontsize = 12
        ax.set_xlabel('Blade thickness', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel('$u$ parameter', fontsize=fontsize, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.axis('off')
        u_plot = np.linspace(0.00, 1.00, 500)

        # Draw thickness distribution
        line, = ax.plot(u_plot, np.real(self.get_upper_thickness_distribution_values(u_plot)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        line, = ax.plot(u_plot, -np.real(self.get_lower_thickness_distribution_values(u_plot)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        # Draw the control polyline
        n = np.shape(self.upper_thickness_distribution.P)[1]
        u_thickness = np.linspace(0, 1, n)
        line, = ax.plot(u_thickness, np.real(self.upper_thickness_distribution.P[0, :]))
        line.set_linewidth(1.25)
        line.set_linestyle("-.")
        line.set_color("r")
        line.set_marker("o")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        n = np.shape(self.lower_thickness_distribution.P)[1]
        u_thickness = np.linspace(0, 1, n)
        line, = ax.plot(u_thickness, -np.real(self.lower_thickness_distribution.P[0, :]))
        line.set_linewidth(1.25)
        line.set_linestyle("-.")
        line.set_color("r")
        line.set_marker("o")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        # Draw thickness distribution at the sample points
        line, = ax.plot(self.u_sample_points, np.real(self.get_upper_thickness_distribution_values(self.u_sample_points)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle(" ")
        line.set_color("k")
        line.set_marker("o")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("b")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        line, = ax.plot(self.u_sample_points, -np.real(self.get_lower_thickness_distribution_values(self.u_sample_points)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle(" ")
        line.set_color("k")
        line.set_marker("o")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("b")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        # # Set the aspect ratio of the data
        # ax.set_aspect(1.0)

        # Set the aspect ratio of the figure
        ratio = 1.00
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        return fig, ax






    # # ---------------------------------------------------------------------------------------------------------------- #
    # # Cascade opening functions (deprecated code)
    # # ---------------------------------------------------------------------------------------------------------------- #
    # def get_cascade_opening(self):
    #
    #     """ Get the cascade opening
    #
    #     Find the cascade opening solving the system of equations:
    #
    #       1) error_x(u,v) = 0
    #       2) error_y(u,v) = 0
    #
    #       where:
    #
    #       (error_x, error_y) = [(x_out, y_out) - v * (n_x, n_y)] - [(x_upper(u), y_upper(u))]
    #
    #     The system of equations is solved with the function fsolve from the scipy library
    #
    #     """
    #
    #     # Initial guess for (u,v)
    #     u0 = 0.7
    #     v0 = 0.1 * self.chord
    #
    #     # Solve the system of equations
    #     solution = fsolve(func=self.get_opening_error,
    #                       x0=np.asarray([u0, v0]),
    #                       args=(),
    #                       fprime=None,
    #                       xtol=1e-8,
    #                       maxfev=1000)
    #
    #     return self.opening
    #
    #
    # def get_opening_error(self, x):
    #
    #     """ Get the error of the current iteration """
    #
    #     # Unpack the independent variables
    #     u = x[0]
    #     v = x[1]
    #
    #     # Get the camberline coordinates and slope at the trailing edge
    #     r_end = self.get_camberline_coordinates(u=1.00)
    #
    #     # Get the vector normal at the trailing edge
    #     normal_end = np.squeeze(self.get_camberline_normal(u=1.00))
    #
    #     # Compute the equation system residual
    #     r1 = r_end - normal_end * v
    #     r2 = self.get_upper_surface_coordinates(u) - np.asarray([0, self.spacing])
    #     error = (r1 - r2).flatten()
    #
    #     # Get opening
    #     self.opening = np.sum((r_end - r1) ** 2) ** (1 / 2)
    #     self.opening2spacing = self.opening / self.spacing
    #     self.normal_end = normal_end
    #     self.AR = self.spacing * np.cos(self.theta_in) / self.opening
    #
    #     return error