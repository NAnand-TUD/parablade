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
from mpl_toolkits.mplot3d import Axes3D


# -------------------------------------------------------------------------------------------------------------------- #
# Define the B-spline curve class
# -------------------------------------------------------------------------------------------------------------------- #
class BSplineCurve:

    """ Create a B-Spline curve object

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

        p : int
            Degree of the B-Spline curve

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped B-Spline

        References
        ----------
        The NURBS book. Chapters 2 and 3
        L. Piegl and W. Tiller
        Springer, second edition

    """

    def __init__(self, P, p, U):

        # Declare input variables as instance variables
        self.P = P
        self.p = p
        self.U = U

        # Declare additional variables as instance variables
        self.u = None
        self.Nu = None
        self.C = None
        self.dC = {}
        self.dP = {}
        self.tangent = None
        self.normal = None
        self.binormal = None
        self.curvature = None



    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute B-Spline curve coordinates
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_BSplineCurve_value(self, u):

        """ Evaluate the coordinates of the B-Spline curve corresponding to the u-parametrization

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline curve

        Returns
        -------
        C : ndarray with shape (ndim, N)
            Array containing B-Spline curve coordinates
            The first dimension of ´C´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´C´ spans the ´u´ parametrization sample points

        """

        # Get the number of sampling points
        if np.isscalar(u) == 1:
            self.Nu = 1
        else:
            self.Nu = u.size

        # Evaluate the spline curve for input u-parametrization
        self.u = u
        self.C = self.evaluate_BSplineCurve(self.P, self.p, self.U, u)

        return self.C


    def evaluate_BSplineCurve(self, P, p, U, u):

        """ Evaluate the coordinates of the B-Spline curve corresponding to the u-parametrization

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

        p : int
            Degree of the B-Spline curve

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped B-Spline

        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline curve

        Returns
        -------
        C : ndarray with shape (ndim, N)
            Array containing B-Spline curve coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the u-parametrization sample points

        """

         # Check the shape of the input parameters
        if P.ndim > 2:
            raise Exception('P must be an array of shape (ndim, n+1)')

        if not np.isscalar(p):
            raise Exception('p must be an scalar')

        if U.ndim > 1:
            raise Exception('U must be an array of shape (r+1=n+p+2,)')

        if np.isscalar(u):
            pass
        elif u.ndim > 1:
            raise Exception('u must be a scalar or an array of shape (N,)')

        # Maximum index of the control points (counting from zero)
        n = np.shape(P)[1] - 1

        # Compute the B-Spline basis polynomials
        N_basis = self.get_basis_polynomials(n, p, U, u)

        # Compute the coordinates of the B-Spline
        # The summations over n is performed exploiting matrix multiplication (vectorized code)
        C = np.dot(P, N_basis)

        return C


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute basis polynomials
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_basis_polynomials(self, n, p, U, u):

        """ Evaluate the n-th B-Spline basis polynomials of degree ´p´ for the input parameter ´u´

        Parameters
        ----------
        n : int
            Number of basis polynomials

        p : int
            Degree of the basis polynomials

        U : ndarray with shape (r+1=n+p+2,)
            Knot vector of the basis polynomials
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped spline

        u : scalar or ndarray with shape (Nu,)
            Scalar or array containing the u-parameter used to evaluate the basis polynomials

        Returns
        -------
        N : ndarray with shape (n, Nu)
            Array containing the basis spline polynomials of order ´p´ evaluated at ´u´
            The first dimension of ´S´ spans the n-th polynomials
            The second dimension of ´S´ spans the u-direction

        """

        # Check the order and the number of basis basis polynomials
        assert (p <= n), 'The order of the basis polynomials cannot be larger than the number of control points minus one'

        # Preliminary computations
        Nu = self.Nu                        # Number of points where the spline curve is evaluated
        m = n + p + 1                       # Number of basis polynomials in the current step of the recursion
        N = np.zeros((m, np.size(u)))       # Initialize the array of basis polynomials

        # First step of the recursion formula (p = 0)
        # The term that is not part of the definition is added to have U[i] <= u <= U[i + 1] when i=n and u=1
        # This problem is mentioned in the NURBS book section 2.5 - Computational Algorithms
        for i in range(m):
            # N[i, :] = 0.0 + 1.0 * (u >= U[i]) * (u <= U[i + 1])
            N[i, :] = 0.0 + 1.0 * (u >= U[i]) * (u < U[i + 1]) + 1.00 * (np.logical_and(u == 1, i == n))

        # Second and next steps of the recursion formula (p = 1, 2, ...)
        p = 1
        while m > n + 1:

            # Update the number of polynomials
            m = m - 1

            # Initialize a temporal variable for the basis polynomial
            # Set data type to complex to avoid problems when computing derivatives (complex step)
            N_bis = np.zeros((m, Nu), dtype=complex)
            # N_bis = np.zeros((m, Nu))

            for i in range(0, m):

                # Compute first factor (avoid division by zero by convention)
                if (U[i + p] - U[i]) == 0:
                    n1 = np.zeros(Nu)
                else:
                    n1 = (u - U[i]) / (U[i + p] - U[i]) * N[i, :]

                # Compute second factor (avoid division by zero by convention)
                if (U[i + p + 1] - U[i + 1]) == 0:
                    n2 = np.zeros(Nu)
                else:
                    n2 = (U[i + p + 1] - u) / (U[i + p + 1] - U[i + 1]) * N[i + 1, :]

                # Compute basis polynomial (recursion formula)
                N_bis[i, :] = n1 + n2

            # Update the order of the polynomials
            p = p + 1

            # Update the array of basis polynomials
            N = N_bis

        return N


    # ---------------------------------------------------------------------------------------------------------------- #
    # Evaluate the derivatives of the B-Spline
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_BSplineCurve_derivative(self, u, derivative_order=1):

        """ Create and evaluate the derivative of the B-Spline curve

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline curve

        derivative_order : int
            Order of the B-Spline curve derivative.
            Set order=0 to get the original B-Spline curve

        Returns
        -------
        C : dictionary of ndarray with shape (ndim, N)
            Dictionary of arrays containing B-Spline curve derivatives
            The keyword is the order of the derivative (0-derivative, 1-derivative, 2-derivative, etc)
            The first dimension of the arrays spans the ´(x,y,z)´ coordinates
            The second dimension of arrays spans the u-parametrization sample points

        """

        # Get the number of sampling points
        if np.isscalar(u) == 1:
            self.Nu = 1
        else:
            self.Nu = u.size

        # Rename the parameters that define the original B-spline
        p = self.p
        P = self.P
        U = self.U

        # Get the parameters that define the order-th derivative of the B-Spline
        # (the order of the implementation matters, be careful when making changes)
        for i in range(derivative_order):
            n = np.shape(P)[1] - 1
            delta_U = U[p+1:n+p+1] - U[1:n+1]
            P = p / delta_U * (P[:, 1:] - P[:, 0:-1])
            p = p - 1
            U = U[1:-1]

        # Evaluate the derivative of the spline curve for input u-parametrization
        dC = self.evaluate_BSplineCurve(P, p, U, u)

        # Store the current derivative in a dictionary
        self.dC[str(derivative_order)] = dC
        self.dP[str(derivative_order)] = P

        return dC


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the curvature of the B-Spline curve
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_BSplineCurve_curvature(self, u):

        # Compute the partial derivatives
        dC = self.get_BSplineCurve_derivative(u, derivative_order=1)
        ddC = self.get_BSplineCurve_derivative(u, derivative_order=2)

        # Evaluate curvature
        numerator = np.sum(np.cross(ddC, dC, axisa=0, axisb=0, axisc=0)**2, axis=0)**(1/2)
        denominator = (np.sum(dC ** 2, axis=0)) ** (3 / 2)
        self.curvature = numerator/denominator

        return self.curvature


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the B-Spline curve Frenet-Serret unitary vectors
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_BSplineCurve_tangent(self, u):

        # Compute the curve derivatives
        dC = self.get_BSplineCurve_derivative(u, derivative_order=1)

        # Evaluate curvature
        numerator = dC
        denominator = (np.sum(numerator ** 2, axis=0)) ** (1 / 2)
        self.tangent = numerator/denominator

        return self.tangent


    def get_BSplineCurve_normal(self, u):

        # Compute the partial derivatives
        dC = self.get_BSplineCurve_derivative(u, derivative_order=1)
        ddC = self.get_BSplineCurve_derivative(u, derivative_order=2)

        # Evaluate curvature
        numerator = np.cross(dC, np.cross(ddC, dC, axisa=0, axisb=0, axisc=0), axisa=0, axisb=0, axisc=0)
        denominator = (np.sum(numerator ** 2, axis=0)) ** (1 / 2)
        self.normal = numerator/denominator

        return self.normal


    def get_BSplineCurve_binormal(self, u):

        # Compute the partial derivatives
        dC = self.get_BSplineCurve_derivative(u, derivative_order=1)
        ddC = self.get_BSplineCurve_derivative(u, derivative_order=2)

        # Evaluate curvature
        numerator = np.cross(dC, ddC, axisa=0, axisb=0, axisc=0)
        denominator = (np.sum(numerator ** 2, axis=0)) ** (1 / 2)
        self.normal = numerator/denominator

        return self.normal


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting function
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_BSplineCurve(self, fig=None, ax=None, curve='yes', control_points='yes', derivative_order=0,
                          color='k', style = '-'):

        # Compute coordinates
        u = np.linspace(0, 1, 500)
        C = np.real(self.get_BSplineCurve_derivative(u, derivative_order))
        P = np.real(self.dP[str(derivative_order)])

        # Number of dimensions of the problem
        ndim = np.shape(self.P)[0]

        if ndim == 1:

            # Create the figure
            if fig is None:
                fig = plt.figure(figsize=(6, 5))
                ax = fig.add_subplot(111)
            ax.set_xlabel('$u$ parameter', fontsize=12, color='k', labelpad=12)
            ax.set_ylabel('B-Spline value', fontsize=12, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            # ax.set_xticks([])
            # ax.set_yticks([])
            # ax.axis('off')

            line, = ax.plot(u, C[0, :])
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
            # ax.set_aspect(1.0)

            # # Set the aspect ratio of the figure
            # ratio = 1.00
            # x1, x2 = ax.get_xlim()
            # y1, y2 = ax.get_ylim()
            # ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

            # Adjust pad
            plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)


        elif ndim == 2:

            # Create the figure
            if fig is None:
                fig = plt.figure(figsize=(6, 5))
                ax = fig.add_subplot(111)

            ax.set_xlabel('$x$ axis', fontsize=12, color='k', labelpad=12)
            ax.set_ylabel('$y$ axis', fontsize=12, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            # ax.set_xticks([])
            # ax.set_yticks([])
            # ax.axis('off')

            if curve == 'yes':
                line, = ax.plot(C[0, :], C[1, :])
                line.set_linewidth(1.25)
                line.set_linestyle(style)
                line.set_color(color)
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if control_points == 'yes':
                line, = ax.plot(P[0, :], P[1, :])
                line.set_linewidth(1.00)
                line.set_linestyle("-.")
                line.set_color("r")
                line.set_marker("o")
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


        elif ndim == 3:

            # Create the figure
            if fig is None:
                fig = mpl.pyplot.figure(figsize=(6, 5))
                ax = fig.add_subplot(111, projection='3d')
            ax.view_init(azim=150, elev=20)
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
            ax.set_xlabel('$x$ axis', fontsize=12, color='k', labelpad=12)
            ax.set_ylabel('$y$ axis', fontsize=12, color='k', labelpad=12)
            ax.set_zlabel('$z$ axis', fontsize=12, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(12)
            ax.xaxis.set_rotate_label(False)
            ax.yaxis.set_rotate_label(False)
            ax.zaxis.set_rotate_label(False)
            # ax.set_xticks([])
            # ax.set_yticks([])
            # ax.set_zticks([])
            # ax.axis('off')

            if curve == 'yes':
                line, = ax.plot(C[0, :], C[1, :], C[2, :])
                line.set_linewidth(1.25)
                line.set_linestyle(style)
                line.set_color(color)
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if control_points == 'yes':
                line, = ax.plot(P[0, :], P[1, :], P[2, :])
                line.set_linewidth(1.00)
                line.set_linestyle("-.")
                line.set_color("r")
                line.set_marker("o")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("r")
                line.set_markerfacecolor("w")
                line.set_label(' ')

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

            # Adjust pad
            plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)


        else: raise Exception('The number of dimensions must be 1, 2 or 3')

        return fig, ax


    def plot_BSplineCurve_curvature(self, fig=None, ax=None):

        # Create the figure
        if fig is None:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)
        ax.set_xlabel('$u$ parameter', fontsize=12, color='k', labelpad=12)
        ax.set_ylabel('Curvature', fontsize=12, color='k', labelpad=12)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(12)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(12)
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.axis('off')

        # Plot the curvature distribution
        u = np.linspace(0, 1, 500)
        curvature = np.real(self.get_BSplineCurve_curvature(u))
        line, = ax.plot(u, curvature)
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
        # ax.set_aspect(1.0)

        # Set the aspect ratio of the figure
        ratio = 1.00
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.set_aspect(np.abs((x2-x1)/(y2-y1))*ratio)

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        return fig, ax