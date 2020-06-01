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
import gmsh_api
import gmsh_api.gmsh as gmsh


#----------------------------------------------------------------------------------------------------------------------#
# Import user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')
from bspline_curve import BSplineCurve


# -------------------------------------------------------------------------------------------------------------------- #
# Define the 2D blade class
# -------------------------------------------------------------------------------------------------------------------- #
class BladeGeo:

    """ Create a 2D blade section object

    The parametrization is based on two arcs (upper surface and lower surface)
    The arcs are constructed imposing a thickness distribution normal to the camber line
    The thickness distributions of the upper and lower surfaces are independent
    The connection between upper and lower surfaces at the leading/trailing edges is G2 continuous

    Parameters
    ----------
    IN_file : string
    Path to the configuration file with the design variable values

        - Axial chord
        - Spacing
        - Stagger angle
        - Inlet metal angle
        - Outlet metal angle
        - Leading edge radius
        - Trailing edge radius
        - Camber line leading edge tangent proportion
        - Camber line trailing edge tangent proportion
        - Upper surface thickness distribution control points
        - Lower surface thickness distribution control points

    """

    # Declare additional variables as instance variables
    u =                             None
    N_points =                      None
    u_sample_points =               None
    camberline_sample_points =      None
    camberline =                    None
    upper_thickness_distribution =  None
    lower_thickness_distribution =  None
    upper_surface_BSpline =         None
    lower_surface_BSpline =         None
    section_coordinates =           None
    upper_periodic_boundary =       None
    lower_periodic_boundary =       None
    inflow_boundary =               None
    outflow_boundary =              None


    def __init__(self, IN):

        # Load an re-scale the design variables
        self.IN = IN
        self.chord_axial = IN["chord_axial"][0]
        self.scale_DVs()

        # Load the blade section variables
        self.spacing = IN["spacing"]
        self.stagger = IN["stagger"]
        self.theta_in = IN["theta_in"]
        self.theta_out = IN["theta_out"]
        self.radius_in = IN["radius_in"]
        self.radius_out = IN["radius_out"]
        self.dist_in = IN["dist_in"]
        self.dist_out = IN["dist_out"]

        self.thickness_upper = [self.radius_in,
                                IN['thickness_upper_1'],
                                IN['thickness_upper_2'],
                                IN['thickness_upper_3'],
                                IN['thickness_upper_4'],
                                IN['thickness_upper_5'],
                                IN['thickness_upper_6'],
                                self.radius_out]

        self.thickness_lower = [self.radius_in,
                                IN['thickness_lower_1'],
                                IN['thickness_lower_2'],
                                IN['thickness_lower_3'],
                                IN['thickness_lower_4'],
                                IN['thickness_lower_5'],
                                IN['thickness_lower_6'],
                                self.radius_out]

        # Set the leading edge at (x,y)=(0,0)
        self.y_in = 0
        self.x_in = 0

        # Set the chord from the axial chord and stagger angle
        self.chord = self.chord_axial / np.cos(self.stagger)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Generate the blade geometry
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_geometry(self):

        # Create camber line
        self.make_camberline()

        # Create thickness distribution
        self.make_thickness_distribution(side_name='upper')
        self.make_thickness_distribution(side_name='lower')

        # Get the sampling points to impose the thickness distribution
        self.N_sample = 10
        self.u_sample_points = self.get_sampling_points(N_samples=self.N_sample, sampling_mode='uniform')

        # Create blade surfaces
        self.make_upper_side()
        self.make_lower_side()

        # Create the boundaries of the flow domain
        self.make_periodic_boundary(boundary_side='upper')
        self.make_periodic_boundary(boundary_side='lower')
        # self.make_periodic_boundary_bis(boundary_side='upper')
        # self.make_periodic_boundary_bis(boundary_side='lower')
        self.make_inflow_boundary()
        self.make_outflow_boundary()


    # ---------------------------------------------------------------------------------------------------------------- #
    # Scale the design variables
    # ---------------------------------------------------------------------------------------------------------------- #
    def scale_DVs(self):

        # Convert each singleton lists into scalars
        # Convert all angles from degrees to radians
        # Scale all lengths using the axial chord
        for variable_name in self.IN:

            if variable_name in ['theta_in', 'theta_out', 'stagger']:
                self.IN[variable_name] = self.IN[variable_name][0]*np.pi/180

            elif variable_name in ['spacing', 'radius_in', 'radius_out', 'dist_in', 'dist_out']:
                self.IN[variable_name] = self.IN[variable_name][0] * self.chord_axial

            elif variable_name in ['thickness_upper_1', 'thickness_upper_2', 'thickness_upper_3',
                                   'thickness_upper_4', 'thickness_upper_5', 'thickness_upper_6',
                                   'thickness_lower_1', 'thickness_lower_2', 'thickness_lower_3',
                                   'thickness_lower_4', 'thickness_lower_5', 'thickness_lower_6']:
                self.IN[variable_name] = self.IN[variable_name][0] * self.chord_axial


    # ---------------------------------------------------------------------------------------------------------------- #
    # Get a distribution of sample points
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_sampling_points(self, N_samples, sampling_mode='uniform'):

        """ Create a 1D distribution of sample points

        Parameters
        ----------
        N_samples : scalar
            Number of sampling points

        sampling_mode : string
            Type of clustering used for the sampling


        Returns
        -------
        sample_points : ndarray with shape (N_sample,)
            Array containing the sample point values

        Notes
        -----
        The following options are available for the variable ´sampling_mode´:

            - Uniform spacing of the sampling points
            - Sigmoid cosine distribution that clusters the sample points towards zero and one
            - Sigmoid parametric distribution that clusters the sample points towards zero and one

        """

        # Default uniform sampling
        offset = 1e-12
        u = np.linspace(0 + offset, 1 - offset, N_samples)

        if sampling_mode == 'uniform':
            sample_points = u

        elif sampling_mode == 'cosine':
            sample_points = np.cos(np.pi / 2 * (1 - u)) ** 2

        elif sampling_mode == 'cluster':
            beta = 3 / 2  # Set beta > 0  | Larger values of beta increase clustering
            sample_points = 1 / (1 + (u / (1 - u)) ** (-beta))

        else:
            raise Exception("Choose a valid option for the sampling distribution: 'uniform', 'cosine', 'cluster'")

        return sample_points


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the camber line
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_camberline(self):

        """ Create the camber line curve (B-Spline of degree 3) """

        # Control point coordintes
        x0 = self.x_in
        y0 = self.y_in
        x1 = self.x_in + self.dist_in * np.cos(self.theta_in)
        y1 = self.y_in + self.dist_in * np.sin(self.theta_in)
        x3 = self.x_in + self.chord * np.cos(self.stagger)
        y3 = self.y_in + self.chord * np.sin(self.stagger)
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

    def get_camberline_normal(self, u):

        """ Evaluate and return the unitary vector normal to the camber line """
        dC = self.camberline.get_BSplineCurve_derivative(u, derivative_order=1)
        tangent = dC / np.sum(dC ** 2, axis=0) ** (1 / 2)
        camberline_normal = np.asarray([-tangent[1, :], tangent[0, :]])

        return camberline_normal


    # ---------------------------------------------------------------------------------------------------------------- #
    # Create the upper and lower thickness distributions
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_thickness_distribution(self, side_name):

        """ Create the thickness distribution (B-Spline of degree 3) """

        # Array of control points
        if side_name == 'upper': P = np.asarray([self.thickness_upper])
        elif side_name == 'lower': P = np.asarray([self.thickness_lower])
        else: raise Exception("Choose a valid side name: 'upper' or 'lower'")

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
        if side_name == 'upper':
            self.upper_thickness_distribution = BSplineCurve(P, p, U)
        elif side_name == 'lower':
            self.lower_thickness_distribution = BSplineCurve(P, p, U)


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

        # Definition of the knot vector (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Get the upper surface set of control points
        camberline_sample = self.camberline.get_BSplineCurve_value(self.u_sample_points)
        normal_sample = +self.get_camberline_normal(self.u_sample_points) # positive sign
        thickness_sample = self.upper_thickness_distribution.get_BSplineCurve_value(self.u_sample_points)
        P = camberline_sample + normal_sample * thickness_sample

        # Collapse the first and last control points into the camber line
        P[:, 0] = camberline_sample[:, 0]
        P[:, -1] = camberline_sample[:, -1]

        # Get additional control points to ensure metal angles an G2 continuity
        P1 = self.get_start_G2_control_point(P, p, U, normal_sample[:, 0], 1 / self.radius_in)
        P2 = self.get_end_G2_control_point(P, p, U, normal_sample[:, -1], 1 / self.radius_out)

        # Concatenate all the control points
        P = np.concatenate((P[:, 0, np.newaxis], P1, P[:, 1:-1], P2, P[:, -1, np.newaxis]), axis=1)

        # # Reverse the order of the control points such that the surface is parametrized counter-clockwise for u in [0,1]
        # P = P[:, ::-1]

        # Get the B-Spline curve object
        self.upper_side = BSplineCurve(P, p, U)


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

        # Get the upper surface array initial set of control points
        camberline_sample = self.camberline.get_BSplineCurve_value(self.u_sample_points)
        normal_sample = -self.get_camberline_normal(self.u_sample_points) # negative sign
        thickness_sample = self.lower_thickness_distribution.get_BSplineCurve_value(self.u_sample_points)
        P = camberline_sample + normal_sample * thickness_sample

        # Collapse the first and last control points into the camber line
        P[:, 0] = camberline_sample[:, 0]
        P[:, -1] = camberline_sample[:, -1]

        # Get additional control points to ensure metal angles an G2 continuity
        P1 = self.get_start_G2_control_point(P, p, U, normal_sample[:, 0], 1 / self.radius_in)
        P2 = self.get_end_G2_control_point(P, p, U, normal_sample[:, -1], 1 / self.radius_out)

        # Concatenate all the control points
        P = np.concatenate((P[:, 0, np.newaxis], P1, P[:, 1:-1], P2, P[:, -1, np.newaxis]), axis=1)

        # Get the B-Spline curve object
        self.lower_side = BSplineCurve(P, p, U)


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
        dist_02 = np.sum((P2 - P0) ** 2) ** (1 / 2)
        alpha = np.arccos(np.dot(P2 - P0, normal) / dist_02)
        dist_01 = np.sqrt(b_2 / a_1 ** 2 / curvature * dist_02 * np.sin(alpha))
        P1 = P0 + normal * dist_01
        P1 = P1[:, np.newaxis]

        return P1


    # ---------------------------------------------------------------------------------------------------------------- #
    # Make the periodic boundaries
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_periodic_boundary(self, boundary_side):

        # Define auxiliary parameters
        L = 0.25 * self.chord_axial
        x_in = self.x_in
        y_in = self.y_in
        x_out = self.x_in + self.chord * np.cos(self.stagger)
        y_out = self.y_in + self.chord * np.sin(self.stagger)

        # Inflow region
        x0 = x_in - 2*L
        y0 = y_in - L*np.tan(self.theta_in)
        x1 = x_in - L
        y1 = y_in - L*np.tan(self.theta_in)
        P_a = np.asarray([[x0, y0], [x1, y1]]).transpose()

        # Blade region (sample the camber line)
        u = np.linspace(0, 1, 10)
        P_b = self.camberline.get_BSplineCurve_value(u)

        # Outflow region
        x0 = x_out + L
        y0 = y_out + L*np.tan(self.theta_out)
        x1 = x_out + 2*L
        y1 = y_out + L*np.tan(self.theta_out)
        P_c = np.asarray([[x0, y0], [x1, y1]]).transpose()

        # Combine the arrays of control points
        P = np.concatenate((P_a, P_b, P_c), axis=1)

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

        # Make the periodic boundary
        if boundary_side == 'upper':
            P[1, :] = P[1, :] + self.spacing/2
            self.upper_periodic_boundary = BSplineCurve(P, p, U)
        elif boundary_side == 'lower':
            P[1, :] = P[1, :] - self.spacing / 2
            self.lower_periodic_boundary = BSplineCurve(P, p, U)
        else: raise Exception("Choose a valid `type` option: 'upper' or 'lower'")


    # Alternative with non-axial inlet/outlet
    def make_periodic_boundary_bis(self, boundary_side):

        # Alternative with npn-zero inlet/outlet slope
        # Define auxiliary parameters
        L = 0.25 * self.chord
        x_in = self.x_in
        y_in = self.y_in
        x_out = self.x_in + self.chord * np.cos(self.stagger)
        y_out = self.y_in + self.chord * np.sin(self.stagger)

        # Inflow region
        x0 = x_in - 2.5*L*np.cos(self.theta_in)
        y0 = y_in - 2.5*L*np.sin(self.theta_in)
        P_a = np.asarray([[x0, y0]]).transpose()

        # Blade region (sample the camber line)
        u = np.linspace(0, 1, 10)
        P_b = self.camberline.get_BSplineCurve_value(u)

        # Outflow region
        x0 = x_out + 2.5*L*np.cos(self.theta_out)
        y0 = y_out + 2.5*L*np.sin(self.theta_out)
        P_c = np.asarray([[x0, y0]]).transpose()

        # Combine the arrays of control points
        P = np.concatenate((P_a, P_b, P_c), axis=1)

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

        # Make the periodic boundary
        if boundary_side == 'upper':
            P[1, :] = P[1, :] + self.spacing/2
            self.upper_periodic_boundary = BSplineCurve(P, p, U)
        elif boundary_side == 'lower':
            P[1, :] = P[1, :] - self.spacing / 2
            self.lower_periodic_boundary = BSplineCurve(P, p, U)
        else: raise Exception("Choose a valid `type` option: 'upper' or 'lower'")



    # ---------------------------------------------------------------------------------------------------------------- #
    # Make the inflow boundary
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_inflow_boundary(self):

        # Inflow region
        x0, y0 = self.lower_periodic_boundary.P[:, 0]
        x1, y1 = self.upper_periodic_boundary.P[:, 0]
        P = np.asarray([[x0, y0], [x1, y1]]).transpose()

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 1

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Make the inflow boundary
        self.inflow_boundary = BSplineCurve(P, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Make the outflow boundary
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_outflow_boundary(self):

        # Inflow region
        x0, y0 = self.lower_periodic_boundary.P[:, -1]
        x1, y1 = self.upper_periodic_boundary.P[:, -1]
        P = np.asarray([[x0, y0], [x1, y1]]).transpose()

        # Maximum index of the control points (counting from zero)
        nn = np.shape(P)[1]
        n = nn - 1

        # Define the order of the basis polynomials
        # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
        # Set p = n (number of control points minus one) to obtain a Bezier
        p = 1

        # Weight of the control points
        # W = np.ones((n + 1,))  # Unitary weight

        # Definition of the knot vectors (clamped spline)
        # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
        U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))

        # Make the inflow boundary
        self.outflow_boundary = BSplineCurve(P, p, U)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting the blade section
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_blade_section(self, fig=None, ax=None,
                           upper_side='yes', upper_side_control_points='no',
                           lower_side='yes', lower_side_control_points='no',
                           camberline='no', camberline_control_points='no', camberline_sample_points='no',
                           leading_edge_radius='no', trailing_edge_radius='no', flow_domain='no'):

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
            camberline_coordinates = np.real(self.camberline.get_BSplineCurve_value(u_plot))
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
            upper_side = np.real(self.upper_side.get_BSplineCurve_value(u_plot))
            line, = ax.plot(upper_side[0, :], upper_side[1, :])
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
            lower_side = np.real(self.lower_side.get_BSplineCurve_value(u_plot))
            line, = ax.plot(lower_side[0, :], lower_side[1, :])
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
            line, = ax.plot(np.real(self.lower_side.P[0, :]), np.real(self.lower_side.P[1, :]))
            line.set_linewidth(0.75)
            line.set_linestyle("-.")
            line.set_color("r")
            line.set_marker("o")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("r")
            line.set_markerfacecolor("w")
            line.set_label('o')

        # Draw the sample points on the camber line
        if camberline_sample_points == 'yes':
            camberline_sample = np.real(self.camberline.get_BSplineCurve_value(self.u_sample_points))
            line, = ax.plot(np.real(camberline_sample[0, :]), np.real(camberline_sample[1, :]))
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
            line, = ax.plot(np.real(self.camberline.P[0, :]), np.real(self.camberline.P[1, :]))
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
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.radius_in * np.cos(self.theta_in)
            y_c = self.y_in + self.radius_in * np.sin(self.theta_in)
            x_circle = x_c + self.radius_in * np.cos(theta)
            y_circle = y_c + self.radius_in * np.sin(theta)
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
            theta = np.linspace(0, 2 * np.pi, 250)
            x_c = self.x_in + self.chord * np.cos(self.stagger) - self.radius_out * np.cos(self.theta_out)
            y_c = self.y_in + self.chord * np.sin(self.stagger) - self.radius_out * np.sin(self.theta_out)
            x_circle = x_c + self.radius_out * np.cos(theta)
            y_circle = y_c + self.radius_out * np.sin(theta)
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

        # Draw the flow domain
        if flow_domain == 'yes':
            u_plot = np.linspace(0, 1, 100)
            upper_periodic = np.real(self.upper_periodic_boundary.get_BSplineCurve_value(u_plot))
            lower_periodic = np.real(self.lower_periodic_boundary.get_BSplineCurve_value(u_plot))
            inflow = np.real(self.inflow_boundary.get_BSplineCurve_value(u_plot))
            outflow = np.real(self.outflow_boundary.get_BSplineCurve_value(u_plot))
            lines = []
            lines.append(ax.plot(upper_periodic[0, :], upper_periodic[1, :])[0])
            lines.append(ax.plot(lower_periodic[0, :], lower_periodic[1, :])[0])
            lines.append(ax.plot(inflow[0, :], inflow[1, :])[0])
            lines.append(ax.plot(outflow[0, :], outflow[1, :])[0])
            for line in lines:
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
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


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting the blade row
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_blade_cascade(self, fig=None, ax=None, blade_section='yes', flow_domain='no'):

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

            # Draw upper and lower sides
            if blade_section == 'yes':

                upper_side = np.real(self.upper_side.get_BSplineCurve_value(u_plot))
                line, = ax.plot(upper_side[0, :], upper_side[1, :] + self.spacing * (k - 1))
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

                lower_side = np.real(self.lower_side.get_BSplineCurve_value(u_plot))
                line, = ax.plot(lower_side[0, :], lower_side[1, :] + self.spacing * (k - 1))
                line.set_linewidth(1.25)
                line.set_linestyle("-")
                line.set_color("k")
                line.set_marker(" ")
                line.set_markersize(3.5)
                line.set_markeredgewidth(1)
                line.set_markeredgecolor("k")
                line.set_markerfacecolor("w")
                line.set_label(' ')

            # Draw the flow domain
            if flow_domain == 'yes':
                u = np.linspace(0, 1, 100)
                upper_periodic = np.real(self.upper_periodic_boundary.get_BSplineCurve_value(u))
                lower_periodic = np.real(self.lower_periodic_boundary.get_BSplineCurve_value(u))
                inflow = np.real(self.inflow_boundary.get_BSplineCurve_value(u))
                outflow = np.real(self.outflow_boundary.get_BSplineCurve_value(u))
                lines = []
                lines.append(ax.plot(upper_periodic[0, :], upper_periodic[1, :] + self.spacing * (k - 1))[0])
                lines.append(ax.plot(lower_periodic[0, :], lower_periodic[1, :] + self.spacing * (k - 1))[0])
                lines.append(ax.plot(inflow[0, :], inflow[1, :] + self.spacing * (k - 1))[0])
                lines.append(ax.plot(outflow[0, :], outflow[1, :] + self.spacing * (k - 1))[0])
                for line in lines:
                    line.set_linewidth(1.25)
                    line.set_linestyle("-")
                    line.set_color("k")
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

        # Hide axes
        plt.axis('off')

        return fig, ax


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting the curvature distribution
    # ---------------------------------------------------------------------------------------------------------------- #
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
        u = np.linspace(0.00, 1.00, 1000)

        # Plot upper side curvature
        line, = ax.plot(u, np.real(self.upper_side.get_BSplineCurve_curvature(u)))
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("b")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label('Lower surface')

        # Plot lower side curvature
        line, = ax.plot(u, np.real(self.lower_side.get_BSplineCurve_curvature(u)))
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("r")
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


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plotting the thickness distribution
    # ---------------------------------------------------------------------------------------------------------------- #
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
        line, = ax.plot(u_plot, np.real(self.upper_thickness_distribution.get_BSplineCurve_value(u_plot)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        line, = ax.plot(u_plot, -np.real(self.lower_thickness_distribution.get_BSplineCurve_value(u_plot)[0, :]))
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
        line, = ax.plot(self.u_sample_points, np.real(self.upper_thickness_distribution.get_BSplineCurve_value(self.u_sample_points)[0, :]))
        line.set_linewidth(0.75)
        line.set_linestyle(" ")
        line.set_color("k")
        line.set_marker("o")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("b")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        line, = ax.plot(self.u_sample_points, -np.real(self.lower_thickness_distribution.get_BSplineCurve_value(self.u_sample_points)[0, :]))
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
    #     # Get the camber line coordinates and slope at the trailing edge
    #     r_end = self.camberline.get_BSplineCurve_value(u=1.00)
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




