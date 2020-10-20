###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

############################### FILE NAME: CAD_functions.py ###################################
# =============================================================================================#
# author: Roberto, Nitish Anand                                                               |
#    :PhD Candidates,                                                                         |
#    :Power and Propulsion, Energy Technology,                                                |
#    :NTNU, TU Delft,                                                                         |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#                                                                                             |
# =============================================================================================#

# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy import optimize

from common import sort_2d_list
from common import plotfun_xy

#----------------------------------------------------------------------------------------------------------------------#
# "Cluster mode" imports
#----------------------------------------------------------------------------------------------------------------------#
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except:
    pass



# -------------------------------------------------------------------------------------------------------------------- #
# Compute the arc length of a parametric curve in n-dimensions
# -------------------------------------------------------------------------------------------------------------------- #
def get_arc_length(C_func, t1, t2):

    """ Compute the arc length of a parametric curve ´C(t) = (x_0(t),..., x_n(t))´ using numerical integration

    Parameters
    ----------
    C_func : function returning ndarray with shape (ndim, N)
        Handle to a function that returns the parametric curve coordinates

    t1 : scalar
        Lower limit of integration for the arclength computation

    t2 : scalar
        Upper limit of integration for the arclength computation

    Returns
    -------
    L : scalar
        Arc length of curve C(t) in the interval [t1, t2]

    """

    # Compute the arc length differential using the central finite differences
    # Be careful with the step-size selected, accurary is not critical, but rounding error must not bloe up
    # It is not possible to use the complex step is the result of the arc-length computation is further differentiated
    def get_arc_legth_differential(t, step=1.e-3):
        # dCdt = np.imag(C_func(t + 1j * step)) / step              # dC/dt = (dx_0/dt, ..., dx_n/dt)
        dCdt = (C_func(t + step) - C_func(t - step))/(2.*step)       # dC/dt = (dx_0/dt, ..., dx_n/dt)
        dLdt = np.sqrt(np.sum(dCdt**2., axis=0))                     # dL/dt = [(dx_0/dt)^2 + ... + (dx_n/dt)^2]^(1/2)
        return dLdt

    # Compute the arc length of C(t) in the interval [t1, t2] by numerical integration
    L = integrate.fixed_quad(get_arc_legth_differential, t1, t2, n=10)[0]

    return L

# -------------------------------------------------------------------------------------------------------------------- #
# Compute the derivative of a 1D curve
# -------------------------------------------------------------------------------------------------------------------- #
def get_curve_derivative(x,y):
    
    """ Compute the derivative of a 1D curve

        Parameters
        ----------
        x: x coordinate
        y: y coordinate
            
        Returns
        ----------
        der_curve : ndarray with shape (ndim,)
                    Array containing the curve derivatives at each point
                    
        Author: Ricardo Puente, 09/2020
                r.puente@imperial.ac.uk

    """
    
    interpol_der = interpolate.CubicSpline(x.real,y).derivative(1)
    der_curve = interpol_der(x)
    
    return der_curve



def get_curve_basis(x,y):
    
    """ Compute the local basis vectors of a 1D curve given as point array

        Parameters
        ----------
        x: x coordinate
        y: y coordinate
            
        Returns
        ----------
        tangent : ndarray with shape (ndim, 2)
                  Array containing the curve tangent vectors at each point
                 
        normal  : ndarray with shape (ndim, 2)
                  Array containing the curve normal vectors at each point
                       
        Author: Ricardo Puente, 09/2020

    """
    der_curve = get_curve_derivative(x,y)
    
    m = len(x)
    
    tangent = np.zeros((m,2),dtype=complex)
    normal  = np.zeros((m,2),dtype=complex)
    
    # Normal computed as cross product of the tangent vector and the outside plane vector
    for i in range(0,m):
        dy = der_curve[i]
        imod = (1.+dy*dy)**(-0.5)
        tangent[i][0] = imod
        tangent[i][1] = dy*imod
        normal[i][0]  = -tangent[i][1]
        normal[i][1]  =  tangent[i][0]
    
    return tangent,normal

def get_curve_basis_spline(x,der_interpolant):
    
    """ Compute the local basis vectors of a 1D curve given as a spline interpolant

        Parameters
        ----------
        x: x coordinate to evaluate the interpolant
        der_interpolant : A cubic spline interpolant derivative object
            
        Returns
        ----------
        tangent : ndarray with shape (2, )
                  Array containing the curve tangent vector
                 
        normal  : ndarray with shape (2, )
                  Array containing the curve normal vector
                       
        Author: Ricardo Puente, 09/2020

    """
    
    tangent = np.zeros((2,),dtype=complex)
    normal  = np.zeros((2,),dtype=complex)
    
    dy = der_interpolant(x)
    
    imod = (1.+dy*dy)**(-0.5)
    tangent[0] = imod
    tangent[1] = dy*imod
    normal[0]  = -tangent[1]
    normal[1]  =  tangent[0]
    
    return tangent,normal

# -------------------------------------------------------------------------------------------------------------------- #
# Compute the bisector curve of 2 others numerically
# -------------------------------------------------------------------------------------------------------------------- #
def get_bisectors(curve_1,curve_2,npoints=50):
    
    """ Compute the bisector curve of 2 other curves

        Parameters
        ----------
        curve_1 : ndarray with shape (2, ndim)
                  Array containing one curve
                  
        curve_1 : ndarray with shape (2,ndim)
                  Array containing the other curve
                  
        npoints : The number of discretization points 
            
        Returns
        ----------
        bisector         : ndarray with shape (2, points_resolved)
                           Array containing the bisector curve
                       
        tangent_points_1 : ndarray with shape ( points_resolved,2) containing the upper curve tangency points with the inscribed circle
                   
        tangent_points_2 : ndarray with shape ( points_resolved,2)  containing the lower curve tangency points with the inscribed circle

         Author: Ricardo Puente, 09/2020
    """

    bisector = []
    tangent_points_1 = []
    tangent_points_2 = []
    
    # Create the interpolant spline object. Use only the real part of x (the thickeness may be complex but is computed as a function of real x coordinates)
    y_1_interpolant = interpolate.CubicSpline(curve_1[0].real,curve_1[1])
    y_2_interpolant = interpolate.CubicSpline(curve_2[0].real,curve_2[1])
    
    
    y_1_interpolant_der = y_1_interpolant.derivative(1)
    y_2_interpolant_der = y_2_interpolant.derivative(1)

    # Create a parametric basis for the curves
    s = np.linspace(0., 1.,npoints)
    
    # Define the function to find the root of
    def equidistance_condition(s,x0,L,xp,yp,np,interpolant,der_interpolant):
        x = x0+L*s
        y = interpolant(x)
        t,n = get_curve_basis_spline(x,der_interpolant)
        dx  = xp-x
        dy  = yp-y
        sn  = np+n
        res = dy*sn[0]-dx*sn[1]     
        return res.real
    
    
    # Data necessary for the curve parametrization
    min_x1 = curve_1[0][0]
    min_x2 = curve_2[0][0]
    
    L1 = curve_1[0][-1]-min_x1
    L2 = curve_2[0][-1]-min_x2
    
    
    ## DBG
    # x1_interpolated = np.zeros(npoints)
    # x2_interpolated = np.zeros(npoints)
    # for i in range(0,npoints):
    #     si = s[i]
    #     x1_interpolated[i] = min_x1+L1*si
    #     x2_interpolated[i] = min_x2+L2*si
    # y1_interpolated = y_1_interpolant(x1_interpolated)
    # y2_interpolated = y_2_interpolant(x2_interpolated) 
    # plotfun_xy(x1_interpolated,y1_interpolated,'y1')
    # plotfun_xy(x2_interpolated,y2_interpolated,'y2')
    ##
    
    
    # # Compute curve lengths        
    # D1 = get_arc_length(y_1_interpolant, min_x1,min_x1+L1)
    # D2 = get_arc_length(y_2_interpolant, min_x2,min_x2+L2)
        
    # fac=D1/D2
    
    for i in range(0,npoints):
        si = s[i]
        # x,y and normals in first curve
        xi = min_x1+L1*si
        yi = y_1_interpolant(xi)
        ti,ni = get_curve_basis_spline(xi,y_1_interpolant_der)
        # Define a single variable lambda
        zero_fun = lambda u: equidistance_condition(u,min_x2,L2,xi,yi,ni,y_2_interpolant,y_2_interpolant_der)
        
        # Find the s in curve 2 where the equidistance condition is fullfilled
        try:
            sol = optimize.root_scalar(zero_fun, bracket=[s[0],s[-1]], method='brentq')
            sj = sol.root
            # x,y and normals in second curve
            xj = min_x2+L2*sj
            yj = y_2_interpolant(xj)        
            tj,nj = get_curve_basis_spline(xj,y_2_interpolant_der)
            # Compute the radius of the circle
                
            sum_nx = ni[0]+nj[0]
            sum_ny = ni[1]+nj[1]
            
            # Prevent zero division with parallel normals
            if np.isclose(sum_nx,0.):
                d = (yi-yj)/sum_ny
            elif np.isclose(sum_ny,0.):
                d = (xi-xj)/sum_nx
            else:
                # These two expressions should give the same d, I take the average for numerical stability
                d1 = (xi-xj)/sum_nx
                d2 = (yi-yj)/sum_ny
                d = 0.5*(d1+d2)
            
            if d < 0.:
                print('Bracketing: ' + str(zero_fun(0.)) +', '+ str(zero_fun(1.)))
                print('Wrong solution to bisection problem!')
            else:
                if d > L1: # Absurd airfoil
                    print('Thickness to chord ratio: ' + str(d/L1))
                    print('Either absurd airfoil or bisection error...')
                else:
                    # Get the bisection point
                    xbi = xi - d * ni[0] 
                    ybi = yi - d * ni[1] 
                    bisector.append([xbi,ybi])
                    # Get the tangent points
                    tangent_points_1.append([xi,yi]) 
                    tangent_points_2.append([xj,yj]) 
        except:
            print('Bracketing: ' + str(zero_fun(0.)) +', '+ str(zero_fun(1.)))
            print('No solution found for bisection problem, not appending...')
        
        
    bisector = np.array(bisector).transpose() 
    tangent_points_1 = np.array(tangent_points_1)
    tangent_points_2 = np.array(tangent_points_2)
    
    return bisector,tangent_points_1,tangent_points_2

# -------------------------------------------------------------------------------------------------------------------- #
# Compute a monotonicity measure of a curve
# -------------------------------------------------------------------------------------------------------------------- #
def get_monotonicity_measure(x,y):
    """ Compute a monotonicity measure of a curve defined by two point arrays
        Ref: "Quantifying non-monotonicity of functions and lack of positivity in signed measures"
              Y. Davydov and R. Zitikis, Modern Stochastics: Theory and Applications
              

        Parameters
        ----------
        x: x coordinate
        y: y coordinate
            
        Returns
        ----------
        LOM : Lack of monotonicity. If LOM <= 0, it is a monotone function
                       
         Author: Ricardo Puente, 09/2020
                r.puente@imperial.ac.uk
    """  

    curve_der_abs = abs(get_curve_derivative(x,y))

    integral = np.trapz(curve_der_abs,x)
    
    LOM = 1.-abs(y[0]-y[-1])/integral
    
    return LOM
 
# -------------------------------------------------------------------------------------------------------------------- #
# Split a closed curve in upper and lower sides
# -------------------------------------------------------------------------------------------------------------------- #   
def split_curve(curve):    
    """ Split a closed curve in upper and lower sides
    
        Parameters
        ----------
        curve_in: Array of size (2,Npoints) with the curve points
            
        Returns
        ----------
        xUp : Upper side
        xDw : Lower side
                       
         Author: Ricardo Puente, 09/2020
                r.puente@imperial.ac.uk
                
    """
    
    idx_mx = np.where(curve[0] == np.amax(curve[0]))[0][0]
    idx_mn = np.where(curve[0] == np.amin(curve[0]))[0][0]
    
    idx1 = min([idx_mx,idx_mn])
    idx2 = max([idx_mx,idx_mn])
    
    x1=[]
    x2=[] 
    for i in range(curve[0].size):
        xi = curve[0][i]
        if i < idx1 or i > idx2:
            x1.append([xi,curve[1][i]])
        else:
            x2.append([xi,curve[1][i]]) 

    # Remove duplicates 
    def remove_duplicates(k):            
        new_k = []
        for elem in k:
            if elem not in new_k:
                new_k.append(elem)            
        return new_k
    
    x1 = remove_duplicates(x1)
    x2 = remove_duplicates(x2)
    
    # Sort the lists according to the x coordinate    
    x1 = sort_2d_list(x1)
    x2 = sort_2d_list(x2) 
    
                    
    xtemp1 = np.array(x1).transpose() 
    xtemp2 = np.array(x2).transpose()
    
    # Check which curve is upper or lower
    max1 = max(xtemp1[1])
    max2 = max(xtemp2[1])
    
    if max1>max2:
        xUp = xtemp1
        xDw = xtemp2
    else:
        xDw = xtemp1
        xUp = xtemp2 
                
    return xUp,xDw


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
        self.dC[str(0) + '-derivative'] = self.C
        self.dP[str(0) + '-derivative'] = self.P

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
            The first dimension of ´C´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´C´ spans the u-parametrization sample points

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
    def get_BSplineCurve_derivative(self, u, order=1):

        """ Create and evaluate the derivative of the B-Spline curve

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline curve

        order : int
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
        for i in range(order):
            n = np.shape(P)[1] - 1
            delta_U = U[p+1:n+p+1] - U[1:n+1]
            P = p / delta_U * (P[:, 1:] - P[:, 0:-1])
            p = p - 1
            U = U[1:-1]

        # Evaluate the derivative of the spline curve for input u-parametrization
        dC = self.evaluate_BSplineCurve(P, p, U, u)

        # Store the current derivative in a dictionary
        self.dC[str(order) + '-derivative'] = dC
        self.dP[str(order) + '-derivative'] = P

        return dC


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_BSplineCurve(self, options, order=0):

        # Choose whether to plot the B-Spline or its derivatives
        C = np.real(self.dC[str(order) + '-derivative'])
        P = np.real(self.dP[str(order) + '-derivative'])

        # Number of dimensions of the problem
        ndim = np.shape(self.P)[0]

        if ndim == 1:
            # Create the figure
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(111)
            fontsize = 12
            ax.set_xlabel('$u$ parameter', fontsize=fontsize, color='k', labelpad=12)
            ax.set_ylabel('B-Spline value', fontsize=fontsize, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
            # ax.set_xticks([])
            # ax.set_yticks([])
            # ax.axis('off')

            line, = ax.plot(self.u, C[0, :])
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


        if ndim == 2:

            # Create the figure
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
            # ax.axis('off')

            if options['line'] == 'yes':
                line, = ax.plot(C[0, :], C[1, :])
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if options['control_points'] == 'yes':
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

            # Prepare the plot
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
            fontsize = 11
            ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=12)
            ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=12)
            ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
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

            if options['line'] == 'yes':
                line, = ax.plot(C[0, :], C[1, :], C[2, :])
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if options['control_points'] == 'yes':
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




# -------------------------------------------------------------------------------------------------------------------- #
# Define the B-spline surface class
# -------------------------------------------------------------------------------------------------------------------- #
class BSplineSurface:

    """ Create a B-Spline surface object

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1, m+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)
            The third dimension of ´P´ spans the v-direction control points (0, 1, ..., m)

        p : int
            Degree of the u-basis polynomials

        q : int
            Degree of the v-basis polynomials

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped spline

        V : ndarray with shape (s+1=m+q+2,)
            The knot vector in the v-direction
            Set the multiplicity of the first and last entries equal to ´q+1´ to obtain a clamped spline

        References
        ----------
        The NURBS book. Chapters 2 and 3
        L. Piegl and W. Tiller
        Springer, second edition

    """

    def __init__(self, P, p, q, U, V):

        # Declare input variables as instance variables
        self.P = P
        self.p = p
        self.q = q
        self.U = U
        self.V = V

        # Declare additional variables as instance variables
        self.u = None
        self.v = None
        self.N = None
        self.S = None


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute B-Spline surface coordinates
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_BSplineSurface_value(self, u, v):

        """ Evaluate the coordinates of the B-Spline surface corresponding to the (u,v) parametrization

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline surface

        v : scalar or ndarray with shape (N,)
            Scalar or array containing the v-parameter used to evaluate the B-Spline surface

        Returns
        -------
        S : ndarray with shape (ndim, N)
            Array containing the B-Spline surface coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the (u,v) parametrization sample points

        """

        # Get the number of u-sampling points
        if np.isscalar(u) == 1:
            Nu = 1
        else:
            Nu = u.size

        # Get the number of v-sampling points
        if np.isscalar(v) == 1:
            Nv = 1
        else:
            Nv = v.size

        # Check that both u and v have the same size
        if Nu == Nv:
            self.N = Nu
        else:
            raise Exception('u and v must have the same size')

        # Evaluate the spline surface for input (u,v) parametrization
        self.S = self.evaluate_BSplineSurface(self.P, self.p, self.q, self.U, self.V, u, v)

        return self.S


    def evaluate_BSplineSurface(self, P, p, q, U, V, u, v):

        """ Evaluate the coordinates of the B-Spline surface corresponding to the (u,v) parametrization

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1, m+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)
            The third dimension of ´P´ spans the v-direction control points (0, 1, ..., m)

        p : int
            Degree of the u-basis polynomials

        q : int
            Degree of the v-basis polynomials

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped spline

        V : ndarray with shape (s+1=m+q+2,)
            The knot vector in the v-direction
            Set the multiplicity of the first and last entries equal to ´q+1´ to obtain a clamped spline

        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the B-Spline surface

        v : scalar or ndarray with shape (N,)
            Scalar or array containing the v-parameter used to evaluate the B-Spline surface

        Returns
        -------
        S : ndarray with shape (ndim, N)
            Array containing the B-Spline surface coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the (u,v) parametrization sample points

        """

        # Check the shape of the input parameters
        if P.ndim > 3:
            raise Exception('P must be an array of shape (ndim, n+1, m+1)')

        if not np.isscalar(p):
            raise Exception('p must be an scalar')

        if not np.isscalar(q):
            raise Exception('q must be an scalar')

        if U.ndim > 1:
            raise Exception('U must be an array of shape (r+1=n+p+2,)')

        if V.ndim > 1:
            raise Exception('V must be an array of shape (s+1=m+q+2,)')

        if np.isscalar(u):
            pass
        elif u.ndim > 1:
            raise Exception('u must be a scalar or an array of shape (N,)')

        if np.isscalar(v):
            pass
        elif u.ndim > 1:
            raise Exception('v must be a scalar or an array of shape (N,)')

        # Shape of the array of control points
        n_dim, nn, mm = np.shape(P)

        # Maximum index of the control points (counting from zero)
        n = nn - 1
        m = mm - 1

        # Compute the B-Spline basis polynomials
        N_basis_u = self.get_basis_polynomials(n, p, U, u)  # shape (n+1, N)
        N_basis_v = self.get_basis_polynomials(m, q, V, v)  # shape (m+1, N)

        # Compute the coordinates of the spline surface (u,v)
        # The implementation of the coordinates computations uses several functions for vectorized code (increase speed)
        A = np.dot(P, N_basis_v)                                        # shape (ndim, n+1, N)
        B = np.repeat(N_basis_u[np.newaxis], repeats=n_dim, axis=0)     # shape (ndim, n+1, N)
        S = np.sum(A*B,axis=1)                                          # shape (ndim, N)

        return S


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
        Nu = self.N                         # Number of points where the spline is evaluated
        m = n + p + 1                       # Number of basis polynomials in the current step of the recursion
        N = np.zeros((m, np.size(u)))       # Initialize the array of basis polynomials

        # First step of the recursion formula (p = 0)
        # The term that is not part of the definition is added to have U[i] <= u <= U[i + 1] when i=n and u=1
        # This problem is mentioned in the NURBS book section 2.5 - Computational Algorithms
        for i in range(m):
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
    # Plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_BSplineSurface(self, options):

        # Prepare the plot
        fig = mpl.pyplot.figure(figsize=(6, 5))
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
        ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=8)
        ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=8)
        ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=8)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        # ax.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
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

        # Get coordinates
        x = np.real(self.S[0, :])
        y = np.real(self.S[1, :])
        z = np.real(self.S[2, :])

        # Get control point coordinates
        Px = np.real(self.P[0, :, :])
        Py = np.real(self.P[1, :, :])
        Pz = np.real(self.P[2, :, :])

        # Plot the spline coordinates as a cloud of points
        if options['point_cloud'] == 'yes':
            ax.plot(x, y, z, linewidth=0, marker='o', markersize='2', markerfacecolor='w')

        # Plot the spline coordinates as a surface
        if options['surface'] == 'yes':
            Nu = options['surface_Nu']
            Nv = options['surface_Nv']
            X = x.reshape(Nv, Nu)
            Y = y.reshape(Nv, Nu)
            Z = z.reshape(Nv, Nu)
            ax.plot_surface(X, Y, Z,
                            color = 'blue',
                            edgecolor = 'black',
                            linewidth = 0.25,
                            alpha = 0.7,
                            shade = True,
                            antialiased = True,
                            zorder = 0)

        # Plot the control mesh coordinates
        if options['control_points'] == 'yes':

            ax.plot_wireframe(Px, Py, Pz,
                              edgecolor = 'red',
                              linewidth = 0.75,
                              linestyles = 'solid',
                              alpha = 1.0,
                              antialiased = True,
                              zorder = 1)

            ax.plot(Px.flatten(), Py.flatten(), Pz.flatten(),
                    linestyle = ' ', linewidth = 1.00, zorder = 4,
                    marker = 'o', markersize = 4, markeredgecolor = 'r', markerfacecolor = 'w', markeredgewidth = 0.75)

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









# -------------------------------------------------------------------------------------------------------------------- #
# Define the NURBS curve class
# -------------------------------------------------------------------------------------------------------------------- #
class NurbsCurve:

    """ Create a NURBS (Non-Uniform Rational Basis Spline) curve object

    Parameters
    ----------
    P : ndarray with shape (ndim, n+1)
        Array containing the coordinates of the control points
        The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
        The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

    W : ndarray with shape (n+1,)
        Array containing the weight of the control points

    p : int
        Degree of the NURBS curve

    U : ndarray with shape (r+1=n+p+2,)
        The knot vector in the u-direction
        Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped NURBS

    References
    ----------
    The NURBS book. Chapter 4
    L. Piegl and W. Tiller
    Springer, second edition

    """

    def __init__(self, P, W, p, U):

        # Declare input variables as instance variables
        self.P = P
        self.W = W
        self.p = p
        self.U = U

        # Declare additional variables as instance variables
        self.u = None
        self.Nu = None
        self.C = None


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute NURBS curve coordinates
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_NurbsCurve_value(self, u):

        """ Evaluate the coordinates of the NURBS curve corresponding to the u-parametrization

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the NURBS curve

        Returns
        -------
        C : ndarray with shape (ndim, N)
            Array containing NURBS curve coordinates
            The first dimension of ´C´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´C´ spans the ´u´ parametrization sample points

        """

        # Get the number of sampling points
        if np.isscalar(u) == 1:
            self.Nu = 1
        else:
            self.Nu = u.size

        # Evaluate the spline curve for input u-parametrization
        self.C = self.evaluate_NurbsCurve(self.P, self.W, self.p, self.U, u)

        return self.C


    def evaluate_NurbsCurve(self, P, W, p, U, u):

        """ Evaluate the coordinates of the NURBS curve corresponding to the u-parametrization

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)

        W : ndarray with shape (n+1,)
            Array containing the weight of the control points

        p : int
            Degree of the NURBS curve

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped NURBS

        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the NURBS curve

        Returns
        -------
        C : ndarray with shape (ndim, N)
            Array containing NURBS curve coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the u-parametrization sample points

        """

        # Check the shape of the input parameters
        if P.ndim > 2:
            raise Exception('P must be an array of shape (ndim, n+1)')

        if W.ndim > 1:
            raise Exception('W must be an array of shape (n+1,)')

        if not np.isscalar(p):
            raise Exception('p must be an scalar')

        if U.ndim > 1:
            raise Exception('U must be an array of shape (r+1=n+p+2,)')

        if np.isscalar(u):
            pass
        elif u.ndim > 1:
            raise Exception('u must be a scalar or an array of shape (N,)')


        # Shape of the array of control points
        n_dim, nn = np.shape(P)

        # Maximum index of the control points (counting from zero)
        n = nn - 1

        # Compute the B-Spline basis polynomials
        N_basis = self.get_basis_polynomials(n, p, U, u)

        # Control points in the expanded space
        # Pw = np.vstack((np.asarray(P*W), W))

        # Map the control points to the expanded space | P = (x*w,y*w,z*w,w)
        w = np.repeat(W[np.newaxis, :], repeats=n_dim, axis=0)
        P_w = np.vstack((np.asarray(P*w), W[np.newaxis, :]))

        # Compute the coordinates of the NURBS curve in the expanded space
        C_w = np.dot(P_w, N_basis)

        # Map the coordinates back to the ordinary space
        C = C_w[0:-1,:]/C_w[-1, :]

        # # Brute force approach to compute the NURBS coordinates without using the expanded space (way slower!)
        # # The implementation of the coordinates computations uses several functions for vectorized code (increase speed)
        # # The details might be easier to understand writing down the NURBS formula and arrays sizes in a piece of paper
        # W = np.repeat(W[:, np.newaxis], repeats=self.Nu, axis=1)
        # R_basis = (N_basis * W)/np.repeat(np.sum(N_basis*W, axis=0)[np.newaxis,:], repeats=n+1, axis=0)
        # C = np.dot(P, R_basis)

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
    # Plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_NurbsCurve(self, options, order=0):

        # Load variables
        C = np.real(self.C)
        P = np.real(self.P)
        ndim = np.shape(self.P)[0]

        if ndim == 2:

            # Create the figure
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
            # ax.axis('off')

            if options['line'] == 'yes':
                line, = ax.plot(C[0, :], C[1, :])
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if options['control_points'] == 'yes':
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

            # Prepare the plot
            fig = plt.figure(figsize=(6, 5))
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
            fontsize = 11
            ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=12)
            ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=12)
            ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=12)
            # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            # ax.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
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

            if options['line'] == 'yes':
                line, = ax.plot(C[0, :], C[1, :], C[2, :])
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            if options['control_points'] == 'yes':
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


# -------------------------------------------------------------------------------------------------------------------- #
# Define the NURBS surface class
# -------------------------------------------------------------------------------------------------------------------- #
class NurbsSurface:

    """ Create a NURBS (Non-Uniform Rational Basis Spline) surface object

    Parameters
    ----------
    P : ndarray with shape (ndim, n+1, m+1)
        Array containing the coordinates of the control points
        The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
        The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)
        The third dimension of ´P´ spans the v-direction control points (0, 1, ..., m)

    W : ndarray with shape (n+1, m+1)
        Array containing the weight of the control points
        The first dimension of ´W´ spans the u-direction control points weights (0, 1, ..., n)
        The second dimension of ´W´ spans the v-direction control points weights (0, 1, ..., m)

    p : int
        Degree of the u-basis polynomials

    q : int
        Degree of the v-basis polynomials

    U : ndarray with shape (r+1=n+p+2,)
        The knot vector in the u-direction
        Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped spline

    V : ndarray with shape (s+1=m+q+2,)
        The knot vector in the v-direction
        Set the multiplicity of the first and last entries equal to ´q+1´ to obtain a clamped spline

    References
    ----------
    The NURBS book. Chapter 4
    L. Piegl and W. Tiller
    Springer, second edition

    """

    def __init__(self, P, W, p, q, U, V):

        # Declare input variables as instance variables
        self.P = P
        self.W = W
        self.p = p
        self.q = q
        self.U = U
        self.V = V

        # Declare additional variables as instance variables
        self.u = None
        self.N = None
        self.v = None
        self.S = None


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute NURBS surface coordinates
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_NurbsSurface_value(self, u, v):

        """ Evaluate the coordinates of the NURBS surface corresponding to the (u,v) parametrization

        Parameters
        ----------
        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the NURBS surface

        v : scalar or ndarray with shape (N,)
            Scalar or array containing the v-parameter used to evaluate the NURBS surface

        Returns
        -------
        S : ndarray with shape (ndim, N)
            Array containing the NURBS surface coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the (u,v) parametrization sample points

        """

        # Get the number of u-sampling points
        if np.isscalar(u) == 1:
            Nu = 1
        else:
            Nu = u.size

        # Get the number of v-sampling points
        if np.isscalar(v) == 1:
            Nv = 1
        else:
            Nv = v.size

        # Check that both u and v have the same size
        if Nu == Nv:
            self.N = Nu
        else:
            raise Exception('u and v must have the same size')

        # Evaluate the spline surface for input (u,v) parametrization
        self.S = self.evaluate_NurbsSurface(self.P, self.W, self.p, self.q, self.U, self.V, u, v)

        return self.S


    def evaluate_NurbsSurface(self, P, W, p, q, U, V, u, v):

        """ Evaluate the coordinates of the NURBS surface corresponding to the (u,v) parametrization

        Parameters
        ----------
        P : ndarray with shape (ndim, n+1, m+1)
            Array containing the coordinates of the control points
            The first dimension of ´P´ spans the coordinates of the control points (1D, 2D, 3D, etc)
            The second dimension of ´P´ spans the u-direction control points (0, 1, ..., n)
            The third dimension of ´P´ spans the v-direction control points (0, 1, ..., m)

        W : ndarray with shape (n+1, m+1)
            Array containing the weight of the control points
            The first dimension of ´W´ spans the u-direction control points weights (0, 1, ..., n)
            The second dimension of ´W´ spans the v-direction control points weights (0, 1, ..., m)

        p : int
            Degree of the u-basis polynomials

        q : int
            Degree of the v-basis polynomials

        U : ndarray with shape (r+1=n+p+2,)
            The knot vector in the u-direction
            Set the multiplicity of the first and last entries equal to ´p+1´ to obtain a clamped spline

        V : ndarray with shape (s+1=m+q+2,)
            The knot vector in the v-direction
            Set the multiplicity of the first and last entries equal to ´q+1´ to obtain a clamped spline

        u : scalar or ndarray with shape (N,)
            Scalar or array containing the u-parameter used to evaluate the NURBS surface

        v : scalar or ndarray with shape (N,)
            Scalar or array containing the v-parameter used to evaluate the NURBS surface

        Returns
        -------
        S : ndarray with shape (ndim, N)
            Array containing the NURBS surface coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the (u,v) parametrization sample points

        """

        # Check the shape of the input parameters
        if P.ndim > 3:
            raise Exception('P must be an array of shape (ndim, n+1, m+1)')

        if W.ndim > 2:
            raise Exception('W must be an array of shape (n+1, m+1)')

        if not np.isscalar(p):
            raise Exception('p must be an scalar')

        if not np.isscalar(q):
            raise Exception('q must be an scalar')

        if U.ndim > 1:
            raise Exception('U must be an array of shape (r+1=n+p+2,)')

        if V.ndim > 1:
            raise Exception('V must be an array of shape (s+1=m+q+2,)')

        if np.isscalar(u):
            pass
        elif u.ndim > 1:
            raise Exception('u must be a scalar or an array of shape (N,)')

        if np.isscalar(v):
            pass
        elif u.ndim > 1:
            raise Exception('v must be a scalar or an array of shape (N,)')

        # Shape of the array of control points
        n_dim, nn, mm = np.shape(P)

        # Maximum index of the control points (counting from zero)
        n = nn - 1
        m = mm - 1

        # Compute the B-Spline basis polynomials
        N_basis_u = self.get_basis_polynomials(n, p, U, u)  # shape (n+1, N)
        N_basis_v = self.get_basis_polynomials(m, q, V, v)  # shape (m+1, N)

        # Map the control points to the expanded space | P = (x*w,y*w,z*w,w)
        w = np.repeat(W[np.newaxis, :, :], repeats=n_dim, axis=0)
        P_w = np.vstack((np.asarray(P*w), W[np.newaxis, :, :]))

        # Compute the coordinates of the NURBS surface in the expanded space
        # The implementation of the coordinates computations uses several functions for vectorized code (increase speed)
        A = np.dot(P_w, N_basis_v)                                      # shape (ndim+1, n+1, N)
        B = np.repeat(N_basis_u[np.newaxis], repeats=n_dim+1, axis=0)   # shape (ndim+1, n+1, N)
        S_w = np.sum(A*B,axis=1)                                        # shape (ndim+1, N)

        # Map the coordinates back to the ordinary space
        S = S_w[0:-1,:]/S_w[-1, :]

        return S


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
        Nu = self.N                         # Number of points where the spline is evaluated
        m = n + p + 1                       # Number of basis polynomials in the current step of the recursion
        N = np.zeros((m, np.size(u)))       # Initialize the array of basis polynomials

        # First step of the recursion formula (p = 0)
        # The term that is not part of the definition is added to have U[i] <= u <= U[i + 1] when i=n and u=1
        # This problem is mentioned in the NURBS book section 2.5 - Computational Algorithms
        for i in range(m):
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
    # Plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_NurbsSurface(self, options):

        """ Plot surface """

        # Prepare the plot
        fig = plt.figure(figsize=(6, 5))
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
        ax.set_xlabel('$x$ axis', fontsize=fontsize, color='k', labelpad=8)
        ax.set_ylabel('$y$ axis', fontsize=fontsize, color='k', labelpad=8)
        ax.set_zlabel('$z$ axis', fontsize=fontsize, color='k', labelpad=8)
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

        # Get coordinates
        x = np.real(self.S[0, :])
        y = np.real(self.S[1, :])
        z = np.real(self.S[2, :])

        # Get control point coordinates
        Px = np.real(self.P[0, :, :])
        Py = np.real(self.P[1, :, :])
        Pz = np.real(self.P[2, :, :])

        # Plot the spline coordinates as a cloud of points
        if options['point_cloud'] == 'yes':
            ax.plot(x, y, z, linewidth=0, marker='o', markersize='2', markerfacecolor='w')

        # Plot the spline coordinates as a surface
        if options['surface'] == 'yes':
            Nu = options['surface_Nu']
            Nv = options['surface_Nv']
            X = x.reshape(Nv, Nu)
            Y = y.reshape(Nv, Nu)
            Z = z.reshape(Nv, Nu)
            ax.plot_surface(X, Y, Z,
                            color = 'blue',
                            edgecolor = 'black',
                            linewidth = 0.25,
                            alpha = 0.7,
                            shade = True,
                            antialiased = True,
                            zorder = 0)

        # Plot the control mesh coordinates
        if options['control_points'] == 'yes':

            ax.plot_wireframe(Px, Py, Pz,
                              edgecolor = 'red',
                              linewidth = 0.75,
                              linestyles = 'solid',
                              alpha = 1.0,
                              antialiased = True,
                              zorder = 1)

            ax.plot(Px.flatten(), Py.flatten(), Pz.flatten(),
                    linestyle = ' ', linewidth = 1.00, zorder = 4,
                    marker = 'o', markersize = 4, markeredgecolor = 'r', markerfacecolor = 'w', markeredgewidth = 0.75)

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


