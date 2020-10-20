###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

################################# FILE NAME: blade_3D.py #####################################
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

#----------------------------------------------------------------------------------------------------------------------#
# Importing general packages
#----------------------------------------------------------------------------------------------------------------------#
import os
import sys
import pdb
import time
import copy
import numpy as np
import scipy.integrate as integrate


#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common import printProgress
from config import ReadUserInput
from CAD_functions import *
from interpolation_functions import *
from blade_2D_connecting_arcs import Blade2DConnectingArcs
from blade_2D_camber_thickness import Blade2DCamberThickness
from common import plotfun 
from common import plotfun_xy
from common import plot3fun

#----------------------------------------------------------------------------------------------------------------------#
# "Cluster mode" imports
#----------------------------------------------------------------------------------------------------------------------#
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except:
    pass


# -------------------------------------------------------------------------------------------------------------------- #
# Define the 3D blade class
# -------------------------------------------------------------------------------------------------------------------- #
class Blade3D:

    """ Create a 3D blade object from configuration file

    Parameters
    ----------
    IN_file : string
        Path to the configuration file that defines the blade geometry

    UV : ndarray with shape (2, N_points)
        Array containing the (u,v) parametrization used to evaluate the blade coordinates
        If UV is not provided Blade3D sets a default (u,v) parametrization


    About the configuration file
    ----------------------------
    The configuration file should specify the following options

    NDIM : Number of dimensions of the problem. Possible values: 2 or 3

    N_BLADES : Number of blades of the cascade. Positive integer

    N_SECTIONS : Number of blade sections used to create the blade (at least 2 sections)
    
    D_STRETCH : Parameter that controls the degree of radial stretching of the sections position using a sigmoid exponential
    
    INTERP_METHOD : Method for 2D surface interpolations. 'bilinear' or 'bicubic' 

    CASCADE_TYPE : Type of blade cascade. Possible values LINEAR or ANNULAR

    PARAMETRIZATION_TYPE : Type of parametrization used to construct the blade sections
                           Set CONNECTING_ARCS to use a 2D parametrization based on connecting arcs
                           Set CAMBER_THICKNESS to use a 2D parametrization based on camberline and thickness

    OPERATION_TYPE : Type of operation executen when calling make_blade().
                     Set GEOMETRY to compute the coordinates of the blade surface
                     Set SENSITIVITY to compute the coordinates, normals and sensitivity of the blade surface

    PLOT_FORMAT : Library used to plot the blade geometry
                  Set MATPLOTLIB to plot the blade using the Matplotlib library
                  Set TECPLOT to plot the blade using the Tecplot 360 Python library


    In addition the configuration file must contain the design variables used to generate the blade geometry

    The shape of the meridional channel is defined giving control points for:
        x_leading, y_leading, z_leading, x_trailing, z_trailing, x_hub, z_hub, x_shroud, z_shroud

    When PARAMETRIZATION_TYPE = CONNECTING_ARCS, the shape of the blade is defined giving control points for:
        stagger, theta_in, theta_out, wedge_in, wedge_out, radius_in, radius_out, dist_1, dist_2, dist_3, dist_4

    When PARAMETRIZATION_TYPE = CAMBER_THICKNESS, the shape of the blade is defined giving control points for:
        stagger, theta_in, theta_out, radius_in, radius_out, dist_in, dist_out, thickness_upper/lower_1...6


    Some hints to write the configuration file
    ------------------------------------------

    It is not necessary to specify string input within quotation marks. Example:

        CASCADE_TYPE = LINEAR       (right)
        CASCADE_TYPE = 'LINEAR'     (wrong)

    It is possible to specify an arbitrary number of control points for any design variable
    To do so, just separate the different values by commas. Example

        stagger = value1, value2, value3

    When contructing the meridional channel, the first and last control points of the hub/shroud are shared with the
    leading and trailing edges. For this reason, these values are not specified in the configuration file and, instead,
    they are computed internally within Blade3D. As a result, if the hub and shroud surfaces are straight lines it is
    not necessary to specify any control point in the configuration file. For example

        x_hub =
        z_hub =
        x_shroud =
        z_shroud =

    would be a valid entries in the configuration file

    In addition, it is not necessary to specify the entry PRESCRIBED_BLADE_FILENAME (see the FitBlade class for info)


    References
    ----------
    TODO add publication information here

    """

    # Array of variable names for the meridional channel
    meridional_channel_names = ['x_leading',
                                'y_leading',
                                'z_leading',
                                'x_trailing',
                                'z_trailing',
                                'x_hub',
                                'z_hub',
                                'x_shroud',
                                'z_shroud']

    # Array of variable names for blade sections (connecting arcs parametrization)
    blade_section_connecting_arcs = ['stagger',
                                   'theta_in',
                                   'theta_out',
                                   'wedge_in',
                                   'wedge_out',
                                   'radius_in',
                                   'radius_out',
                                   'dist_1',
                                   'dist_2',
                                   'dist_3',
                                   'dist_4']

    # Array of variable names for blade sections (camber thickness parametrization)
    blade_section_camber_thickness = ['stagger',
                                       'theta_in',
                                       'theta_out',
                                       'radius_in',
                                       'radius_out',
                                       'dist_in',
                                       'dist_out',
                                       'thickness_upper_1',
                                       'thickness_upper_2',
                                       'thickness_upper_3',
                                       'thickness_upper_4',
                                       'thickness_upper_5',
                                       'thickness_upper_6',
                                       'thickness_lower_1',
                                       'thickness_lower_2',
                                       'thickness_lower_3',
                                       'thickness_lower_4',
                                       'thickness_lower_5',
                                       'thickness_lower_6']

    # Declare additional variables as instance variables
    DVs_names                   = None
    DVs_names_2D                = None
    DVs_names_meridional        = None
    DVs_number                  = None
    DVs_functions               = None
    DVs_control_points          = None
    surface_interpolant         = None
    surface_coordinates         = None
    surface_normals             = None
    surface_sensitivity         = None
    u_hub                       = None
    u_shroud                    = None
    hub_coordinates             = None
    shroud_coordinates          = None
    get_meridional_channel_x    = None
    get_meridional_channel_z    = None
    meanline_length             = None
    Nu                          = None
    Nv                          = None
    section_coordinates         = None

    def __init__(self, IN, UV=None):

        # Declare input variables as instance variables
        self.IN = IN
        self.N_BLADES               = int(IN["N_BLADES"][0])
        self.NDIM                   = int(IN["NDIM"][0])
        self.N_SECTIONS             = int(IN["N_SECTIONS"][0])
        self.CASCADE_TYPE           = IN["CASCADE_TYPE"]
        self.PARAMETRIZATION_TYPE   = IN["PARAMETRIZATION_TYPE"]
        self.OPERATION_TYPE         = IN['OPERATION_TYPE']
        self.PLOT_FORMAT            = IN["PLOT_FORMAT"]
        self.D_STRETCH              = float(IN["D_STRETCH"][0])
         
        #  The interp_method key should be optional
        try:
           self.interp_method = IN["INTERP_METHOD"]
        except:
            self.interp_method = 'bilinear'
            

        # Check the number of sections
        if self.N_SECTIONS < 2: raise Exception('It is necessary to use at least 2 sections to generate the blade')

        # Initialize design variables
        self.initialize_DVs_names()
        self.update_DVs_control_points(IN)
        self.get_DVs_functions()

        # Initialize the (u,v) parametrization
        if UV is None:
            self.initialize_uv_values()
        else:
            self.u = UV[0, :]
            self.v = UV[1, :]
            self.N_points = np.shape(UV)[1]
            
        self.make_surface_interpolant()

    # ---------------------------------------------------------------------------------------------------------------- #
    # Generate the blade geometry and/or sensitities
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_blade(self):

        """ Compute the blade surface coordinates and the sensitivity with respect to the design variables """

        if self.OPERATION_TYPE == 'GEOMETRY':

            # Compute the coordinates of the blade surface
            # self.make_surface_interpolant(interp_method='bilinear')   # TODO bicubic interpolation is not ready
            self.surface_coordinates = self.get_surface_coordinates(self.u, self.v)
            self.make_hub_surface()
            self.make_shroud_surface()

        elif self.OPERATION_TYPE == 'SENSITIVITY':

            # Compute the coordinates of the blade surface
            print("\nStarting coordinates computation...", end='     ')
            # self.make_surface_interpolant(interp_method='bilinear')
            self.surface_coordinates = self.get_surface_coordinates(self.u, self.v)
            self.make_hub_surface()
            self.make_shroud_surface()
            print("Done!")

            # Compute the unitary vectors normal to the blade surface
            print("Starting normals computation...", end='         ')
            self.surface_normals = self.get_surface_normals(self.u, self.v, method='complex_step')
            print("Done!")

            # Compute the derivative of the blade coordinates with respect to the design variables
            print("Starting sensitivity computation...\n")
            self.surface_sensitivity = self.get_surface_sensitivity(self.u, self.v, method='complex_step')
            print('Sensitivity computation is finished!\n')

        else:
            raise Exception('Choose a valid option for OPERATION_TYPE: "GEOMETRY" or "SENSITIVITY"')


    # ---------------------------------------------------------------------------------------------------------------- #
    # Initialize the design variable names
    # ---------------------------------------------------------------------------------------------------------------- #
    def initialize_DVs_names(self):

        """ Initialize the names of the design variables """
        if self.PARAMETRIZATION_TYPE == "CONNECTING_ARCS":
            self.DVs_names = self.meridional_channel_names + self.blade_section_connecting_arcs
            self.DVs_names_2D = self.blade_section_connecting_arcs
        elif self.PARAMETRIZATION_TYPE == "CAMBER_THICKNESS":
            self.DVs_names = self.meridional_channel_names + self.blade_section_camber_thickness
            self.DVs_names_2D = self.blade_section_camber_thickness
        else:
            raise Exception('Choose a valid option for PARAMETRIZATION_TYPE: "CONNECTING_ARCS" or "CAMBER_THICKNESS"')

        # Combine the lists of design variables
        self.DVs_names_meridional = self.meridional_channel_names
        self.DVs_number = len(self.DVs_names)

    # ---------------------------------------------------------------------------------------------------------------- #
    # Get the control points of the design variables
    # ---------------------------------------------------------------------------------------------------------------- #
    def update_DVs_control_points(self, IN):

        """ Create a dictionary with the control points of each design variable """

        # Use a deepcopy to avoid angle conversion problems when calling this method from the BladeFit class
        IN = copy.deepcopy(IN)

        # Assign the values from the IN dictionary to the control_points dictionary
        self.DVs_control_points = {}
        for k in self.DVs_names:
            self.DVs_control_points[k] = IN[k]

        # Adjust the hub and shroud control points so that the shared points match exactly
        self.DVs_control_points['x_hub'] = [IN['x_leading'][0]] + IN['x_hub'] + [IN['x_trailing'][0]]
        self.DVs_control_points['z_hub'] = [IN['z_leading'][0]] + IN['z_hub'] + [IN['z_trailing'][0]]
        self.DVs_control_points['x_shroud'] = [IN['x_leading'][-1]] + IN['x_shroud'] + [IN['x_trailing'][-1]]
        self.DVs_control_points['z_shroud'] = [IN['z_leading'][-1]] + IN['z_shroud'] + [IN['z_trailing'][-1]]

        self.make_surface_interpolant()

    # ---------------------------------------------------------------------------------------------------------------- #
    #  Create the design variable functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_DVs_functions(self):

        """ Create the design variable functions

            Each design variable is described by a law of evolution given by a B-Spline curve of, at most, degree 3
            The design variable functions can be constructed with an arbitrary number of control points (from 1 to N)

        """

        self.DVs_functions = {}
        for k in self.DVs_names:

            # Array of control points
            P = np.array([self.DVs_control_points[k]])

            # Maximum index of the control points (counting from zero)
            nn = np.shape(P)[1]
            n = nn-1

            # Define the order of the basis polynomials
            # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
            # Set p = n (number of control points minus one) to obtain a Bezier
            p = min(n, 3)   # Set at most p = 3 (cubic B-Spline)

            # # Weight of the control points
            # W = np.ones((n+1,))    # Unitary weight for B-Splines

            # Definition of the knot vectors (clamped spline)
            # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
            U = np.concatenate((np.zeros(p), np.linspace(0, 1, n-p+2), np.ones(p)))

            # Get the design variable function object
            self.DVs_functions[k] = BSplineCurve(P, p, U).get_BSplineCurve_value


    # ---------------------------------------------------------------------------------------------------------------- #
    # Initialize the (u,v) parametrization
    # ---------------------------------------------------------------------------------------------------------------- #
    def initialize_uv_values(self, Nu = 500, Nv = 25):

        """ Set a default (u,v) parametrization with the desired number of u- and v-points """

        # print('Using a default (u,v) parametrization...')

        # Set a small offset to ensure that (u,v) are within [0,1] when using finite differences
        h = 1e-5

        # Define a default (u,v) parametrization from a meshgrid
        if self.NDIM == 2:
            self.Nu, self.Nv = Nu, 1
            u = np.linspace(0.00+h, 1.00-h, self.Nu)
            v = 0.50
        elif self.NDIM == 3:
            self.Nu, self.Nv = Nu, Nv
            u = np.linspace(0.00+h, 1.00-h, self.Nu)
            v = np.linspace(0.00+h, 1.00-h, self.Nv)
        else:
            raise Exception('The number of dimensions must be "2" or "3"')

        [u, v] = np.meshgrid(u, v)
        self.u = u.flatten()
        self.v = v.flatten()
        self.N_points = self.Nu*self.Nv


    # ---------------------------------------------------------------------------------------------------------------- #
    # Update the (u,v) parametrization
    # ---------------------------------------------------------------------------------------------------------------- #
    def update_uv_values(self, u, v):

        """ Update the (u,v) parameters with new values given as input """

        Nu = np.size(u)
        Nv = np.size(v)
        if Nu == Nv:
            self.u = u
            self.v = v
            self.N_points = Nu
        else:
            raise Exception('the u- and v-arrays must have the same number of points')


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the coordinates of the blade surface
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_surface_coordinates(self, u, v):

        """ Compute the blade surface coordinates for the current (u,v) parametrization by interpolation

        Parameters
        ----------
        u : ndarray with shape (N,)
            Array containing the u-parameter used to evaluate the blade surface coordinates

        v : ndarray with shape (N,)
            Array containing the v-parameter used to evaluate the blade surface coordinates

        Returns
        -------
        surface_coordinates : ndarray with shape (3, N)
            Array containing the blade surface coordinates
            The first dimension of ´S´ spans the ´(x,y,z)´ coordinates
            The second dimension of ´S´ spans the (u,v) parametrization sample points

        """

        # Check that the query (u,v) parametrization is within the interpolation range
        if np.any(u < 0) or np.any(u > 1) or np.any(v < 0) or np.any(v > 1):
            raise ValueError('Extrapolation of the (u,v) parametrization is not supported')

        # Compute the surface coordinates by interpolation
        surface_coordinates = self.surface_interpolant(u, v)

        return surface_coordinates


    def make_surface_interpolant(self):

        """ Create a surface interpolant using the coordinates of several blade sections

        Set interp_method='bilinear' to create a bilinear interpolator
        Set interp_method='bicubic' to create a bicubic interpolator

        The interpolant created by this method is used by self.get_surface_coordinates() to evaluate the blade surface
        coordinates for any input (u,v) parametrization

        """

        # Update the geometry before creating the interpolant
        self.get_DVs_functions()            # First update the functions that return design variables
        self.make_meridional_channel()      # Then update the functions that return the meridional channel coordinates

        # Compute the coordinates of several blade sections
        # The (u,v) parametrization used to compute the blade sections must be regular (necessary for interpolation)
        num_points_section = 500
        u_interp = np.linspace(0.00, 1., num_points_section)
        v_interp = self.make_radial_distribution()
        S_interp = np.zeros((3, num_points_section, self.N_SECTIONS), dtype=complex)
        self.section_coordinates = np.zeros((self.N_SECTIONS, 2, num_points_section ), dtype=complex)  
        for k in range(self.N_SECTIONS):
            S_interp[..., k],self.section_coordinates[k,...] = self.get_section_coordinates(u_interp, v_interp[k])          
                

        # Create the interpolator objects for the (x,y,z) coordinates
        if self.interp_method == 'bilinear':
            x_function = BilinearInterpolation(u_interp, v_interp, S_interp[0, ...])
            y_function = BilinearInterpolation(u_interp, v_interp, S_interp[1, ...])
            z_function = BilinearInterpolation(u_interp, v_interp, S_interp[2, ...])

        elif self.interp_method == 'bicubic':
            # There seems to be a problem with bicubic interpolation, the output is wiggly
            x_function = BicubicInterpolation(u_interp, v_interp, S_interp[0, ...])
            y_function = BicubicInterpolation(u_interp, v_interp, S_interp[1, ...])
            z_function = BicubicInterpolation(u_interp, v_interp, S_interp[2, ...])

        else:
            raise Exception('Choose a valid interpolation method: "bilinear" or "bicubic"')

        # Create the surface interpolant using a lambda function to combine the (x,y,z) interpolants
        self.surface_interpolant = lambda u, v: np.asarray((x_function(u, v), y_function(u, v), z_function(u, v)))
        
    def make_radial_distribution(self):

        """ Define a radial distribution of sections with an exponential accumulation to the endwalls
        
            Author: Ricardo Puente, 09/2020
                    r.puente@imperial.ac.uk
        """
        
        v_interp = np.linspace(0.,1.,self.N_SECTIONS)
        
        if self.D_STRETCH>0.:
        
            v1 = np.zeros(self.N_SECTIONS)        
        
            x  = np.linspace(-1.,1.,self.N_SECTIONS)
            
            # Compute sigmoid
            e = np.exp(-self.D_STRETCH*x)
            for i in range(self.N_SECTIONS):
                v1[i] = 1./(1.+e[i])
            
            # Normalize for output
            vmn = min(v1)
            iL  = 1./(max(v1)-vmn)
            for i in range(self.N_SECTIONS):
                v_interp[i] = (v1[i]-vmn)*iL
            
        #plotfun_xy(np.linspace(0.,1.,self.N_SECTIONS),v_interp,'radist')
            
        return v_interp


    def make_meridional_channel(self):

        """ Create the functions used to compute the (x,z) cooordinates at any point of the meridional channel

        The coordinates at the interior of the meridional channel are computed by transfinite interpolation from the
        coordinates at the boundaries (leading edge, lower surface, trailing edge and upper surface)

        The coordinates at any interior point can be evaluated as x(u,v) and z(u,v)

        The functions created by this method are used by self.get_section_coordinates() to map the coordinates of
        the 2D unitary blades into the meridional channel

        """

        # Define a function to compute the x-coordinate of the meridional channel as x(u,v)
        self.get_meridional_channel_x = TransfiniteInterpolation(self.DVs_functions["x_leading"],
                                                                 self.DVs_functions["x_hub"],
                                                                 self.DVs_functions["x_trailing"],
                                                                 self.DVs_functions["x_shroud"],
                                                                 self.DVs_control_points["x_leading"][0],
                                                                 self.DVs_control_points["x_trailing"][0],
                                                                 self.DVs_control_points["x_trailing"][-1],
                                                                 self.DVs_control_points["x_leading"][-1])

        # Define a function to compute the z-coordinate of the meridional channel as z(u,v)
        self.get_meridional_channel_z = TransfiniteInterpolation(self.DVs_functions["z_leading"],
                                                                 self.DVs_functions["z_hub"],
                                                                 self.DVs_functions["z_trailing"],
                                                                 self.DVs_functions["z_shroud"],
                                                                 self.DVs_control_points["z_leading"][0],
                                                                 self.DVs_control_points["z_trailing"][0],
                                                                 self.DVs_control_points["z_trailing"][-1],
                                                                 self.DVs_control_points["z_leading"][-1])

        # Compute the arc length of the blade meanline (secondary computation required for the BladeFit class)
        x_func = lambda u: self.get_meridional_channel_x(u, v=0.50)
        z_func = lambda u: self.get_meridional_channel_z(u, v=0.50)
        m_func = lambda u: np.concatenate((x_func(u), z_func(u)), axis=0)
        self.meanline_length = get_arc_length(m_func, 0.0 + 1e-6, 1.0 - 1e-6)


    def get_section_coordinates(self, u_section, v_section):

        """ Compute the coordinates of the current blade section

        Parameters
        ----------
        u_section : ndarray with shape (Nu,)
            Array containing the u-parameter used to evaluate section coordinates

        v_section : scalar
            Scalar containing the v-parameter of the current blade span

        Returns
        -------
        section_coordinates : ndarray with shape (3, Nu)
            Array containing the (x,y,z) coordinates of the current blade section

        plane_section_coordinates : ndarray with shape (2, Nu)
            Array containing the (m',rTh) coordinates of the current blade section

        """
        # Local auxiliary functions
        def transfinite_interpolation(x, v_section):
            # Ensure that the x-coordinates are between zero and one (actually between zero+eps and one-eps)
            uu_section = (x - np.amin(x) + 1e-12)/(np.amax(x) - np.amin(x) + 2e-12)

            # Obtain the x-z coordinates by transfinite interpolation of the meridional channel contour
            xout = self.get_meridional_channel_x(uu_section, v_section)       # x corresponds to the axial direction
            zout = self.get_meridional_channel_z(uu_section, v_section)       # x corresponds to the radial direction

            # Create a single-variable function with the coordinates of the meridional channel
            x_func = lambda u: self.get_meridional_channel_x(u, v_section)
            z_func = lambda u: self.get_meridional_channel_z(u, v_section)
            m_func = lambda u: np.concatenate((x_func(u), z_func(u)), axis=0)

            # Compute the arc length of the meridional channel (streamline)
            arc_length = get_arc_length(m_func, 0.0 + 1e-6, 1.0 - 1e-6)
        
            return xout,zout,arc_length


        # Get the design variables of the current blade section
        section_variables = {}
        for k in self.DVs_names_2D:
            section_variables[k] = self.DVs_functions[k](v_section)

        # Compute the coordinates of a blade section with an unitary meridional chord
        if self.PARAMETRIZATION_TYPE == 'CONNECTING_ARCS':
            section_coordinates = Blade2DConnectingArcs(section_variables).get_section_coordinates(u_section)
        elif self.PARAMETRIZATION_TYPE == 'CAMBER_THICKNESS':
            section_coordinates = Blade2DCamberThickness(section_variables).get_section_coordinates(u_section)
        else:
            raise Exception('Choose a valid option for PARAMETRIZATION_TYPE: "CONNECTING_ARCS" or "CAMBER_THICKNESS"')

        # Rename the section coordinates
        x = section_coordinates[0, :]                   # x corresponds to the meridional direction
        y = section_coordinates[1, :]                   # y corresponds to the tangential direction

        # Ensure that the x-coordinates are between zero and one (actually between zero+eps and one-eps)
        uu_section = (x - np.amin(x) + 1e-12)/(np.amax(x) - np.amin(x) + 2e-12)

        # Obtain the x-z coordinates by transfinite interpolation of the meridional channel contour
        x = self.get_meridional_channel_x(uu_section, v_section)       # x corresponds to the axial direction
        z = self.get_meridional_channel_z(uu_section, v_section)       # x corresponds to the radial direction

        # Create a single-variable function with the coordinates of the meridional channel
        x_func = lambda u: self.get_meridional_channel_x(u, v_section)
        z_func = lambda u: self.get_meridional_channel_z(u, v_section)
        m_func = lambda u: np.concatenate((x_func(u), z_func(u)), axis=0)

        # Compute the arc length of the meridional channel (streamline)
        arc_length = get_arc_length(m_func, 0.0 + 1e-6, 1.0 - 1e-6)

        # Compute the y-coordinates of the current blade section by scaling the unitary blade
        y = self.DVs_functions["y_leading"](v_section) + y * arc_length

        # Piece together the planar x,y coordinates
        plane_section_coordinates = np.concatenate((x,y), axis=0)
    
        # Transform the blade coordinates to cartesian coordinates
        if self.CASCADE_TYPE == "LINEAR":
            X = x
            Y = y
            Z = z
        elif self.CASCADE_TYPE == "ANNULAR":
            X = x
            Y = z*np.sin(y/z)
            Z = z*np.cos(y/z)
        else:
            raise Exception('Choose a valid cascade type: "LINEAR" or "ANNULAR"')

        # Piece together the X-Y-Z coordinates
        section_coordinates = np.concatenate((X, Y, Z), axis=0)

        # DBG
        #plotfun(plane_section_coordinates,'Airfoil')
        ###
        
        return section_coordinates, plane_section_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the unitary vectors normal to the blade surface
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_surface_normals(self, u, v, method='complex_step', step=1e-12):

        """ Compute unitary vectors normal to the blade surface

         Parameters
        ----------
        u : ndarray with shape (N,)
            Array containing the u-parameter used to evaluate the blade surface coordinates

        v : ndarray with shape (N,)
            Array containing the v-parameter used to evaluate the blade surface coordinates

        method : string
            Method used to compute the derivatives of the surface coordinates
            Valid options: 'forward_finite_differences', 'central_finite_differences' or 'complex_step'

        step : scalar
            Step size used to compute the derivatives of the blade surface

        Returns
        -------
        surface_normals : ndarray with shape (3, N)
            Array containing the unitary vectors normal to the blade surface
            The first dimension spans the ´(x,y,z)´ components
            The second dimension spans the (u,v) parametrization sample points

        """

        # Compute a pair of tangent vectors by differentiation of the surface coordinates with respect to (u,v)
        if method == 'forward_finite_differences':
            C = self.get_surface_coordinates(u, v)
            C_u = self.get_surface_coordinates(u + step, v)
            C_v = self.get_surface_coordinates(u, v + step)
            T_u = (C_u - C) / step                              # Forward finite differences in u
            T_v = (C_v - C) / step                              # Forward finite differences in v
            N = np.cross(T_u, T_v, axisa=0, axisb=0, axisc=0)   # Normal vector as cross product of two tangent vectors
            norm = np.sum(N**2, axis=0)**(1/2)                  # 2-norm of the normal vectors
            surface_normals = -N / norm                         # Unitary normal vector (fix sign to point outwards)

        elif method == 'central_finite_differences':
            C_u1 = self.get_surface_coordinates(u - step, v)
            C_u2 = self.get_surface_coordinates(u + step, v)
            C_v1 = self.get_surface_coordinates(u, v - step)
            C_v2 = self.get_surface_coordinates(u, v + step)
            T_u = (C_u2 - C_u1) / (2*step)                      # Central finite differences in u
            T_v = (C_v2 - C_v1) / (2*step)                      # Central finite differences in v
            N = np.cross(T_u, T_v, axisa=0, axisb=0, axisc=0)   # Normal vector as cross product of two tangent vectors
            norm = np.sum(N**2, axis=0)**(1/2)                  # 2-norm of the normal vectors
            surface_normals = -N / norm                         # Unitary normal vector (fix sign to point outwards)

        elif method == 'complex_step':
            C_u = self.get_surface_coordinates(u + step*1j, v)
            C_v = self.get_surface_coordinates(u, v + step*1j)
            T_u = np.imag(C_u) / step                           # Complex step in u
            T_v = np.imag(C_v) / step                           # Complex step in v
            N = np.cross(T_u, T_v, axisa=0, axisb=0, axisc=0)   # Normal vector as cross product of two tangent vectors
            norm = np.sum(N**2, axis=0)**(1/2)                  # 2-norm of the normal vectors
            surface_normals = -N / norm                         # Unitary normal vector (fix sign to point outwards)

        else:
            raise Exception('Choose a valid differentiation method: '
                            '"forward_finite_differences", "central_finite_differences", or "complex_step"')

        return surface_normals


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the derivatives of the surface coordinates with respect to the design variables
    # ---------------------------------------------------------------------------------------------------------------- #
    def get_surface_sensitivity(self, u, v, method='complex_step', step=1e-12, variable='all', display_progress='yes'):

        """ Compute derivative of the blade coordinates with respect to the design variables

         Parameters
        ----------
        u : ndarray with shape (N,)
            Array containing the u-parameter used to evaluate the blade surface coordinates

        v : ndarray with shape (N,)
            Array containing the v-parameter used to evaluate the blade surface coordinates

        method : string
            Method used to compute the derivatives of the surface coordinates
            Valid options: 'forward_finite_differences', 'central_finite_differences' or 'complex_step'

        step : scalar
            Step size used to compute the derivatives of the blade surface

        variable : string or list
            Set variable='all' to compute the sensitivity with respect to all design variables
            Set variable=[name,number] to compute the sensitivity with respect to a single variable (used for debugging)

        display : string
            Set display_progress='yes' to print a bar showing the progress of the sensitivity computation

        Returns
        -------
        surface_sensitivity : dictionary of ndarrays with shape (3, N)
            The keys of the dictionary are the names of the design variables
            The entries of the dictionary are arrays that contain the derivative of the surface coordinates with
            respect to each design variable

        """

        # Store the original set of control points in a deep-copy
        control_points = copy.deepcopy(self.DVs_control_points)

        # Fix the sensitivity overload bug having two dictionaries
        surface_sensitivity = {}
        my_numbers = {}

        # Get the name of the variable to compute derivatives and the number of control points
        if variable == 'all':
            my_keys = self.DVs_names
            for key in my_keys:
                my_numbers[key] = range(len(self.DVs_control_points[key]))
        else:
            my_keys = [variable[0]]
            my_numbers[variable[0]] = [variable[1]]

        # Initialize variables to print progress
        total = len(my_keys)
        count = 0

        # Loop through all the keys and numbers (e.g. key=dist_1, number=0 means dist_1_0)
        for key in my_keys:
            count = count + 1
            for number in my_numbers[key]:

                if method == 'forward_finite_differences':
                    self.DVs_control_points[key][number] = control_points[key][number]
                    self.make_surface_interpolant()
                    C_1 = self.get_surface_coordinates(u, v)
                    self.DVs_control_points[key][number] = control_points[key][number] + step
                    self.make_surface_interpolant()
                    C_2 = self.get_surface_coordinates(u, v)
                    surface_sensitivity[key + '_' + str(number)] = (C_2 - C_1) / step

                elif method == 'central_finite_differences':
                    self.DVs_control_points[key][number] = control_points[key][number] - step
                    self.make_surface_interpolant()
                    C_1 = self.get_surface_coordinates(u, v)
                    self.DVs_control_points[key][number] = control_points[key][number] + step
                    self.make_surface_interpolant()
                    C_2 = self.get_surface_coordinates(u, v)
                    surface_sensitivity[key + '_' + str(number)] = (C_2 - C_1) / (2 * step)

                elif method == 'complex_step':
                    self.DVs_control_points[key][number] = control_points[key][number] + step * 1j
                    self.make_surface_interpolant()
                    C = self.get_surface_coordinates(u, v)
                    surface_sensitivity[key + '_' + str(number)] = np.imag(C) / step

                else:
                    raise Exception('Choose a valid differentiation method: '
                                    '"forward_finite_differences", "central_finite_differences", or "complex_step"')


                # Retrieve the original set of control points
                self.DVs_control_points = copy.deepcopy(control_points)

            # Print progress bar
            if display_progress == 'yes':
                printProgress(count, total)

        # Retrieve the original interpolant
        self.make_surface_interpolant()

        return surface_sensitivity


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the coordinates of the hub surface
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_hub_surface(self):
        self.u_hub = np.linspace(0, 1, 100)
        self.hub_coordinates = self.get_extended_hub_coordinates(self.u_hub)

    def get_extended_hub_coordinates(self, u):

        """ Compute the coordinates of the hub surface in the (x,z) plane

        The hub surface is is extended to the inlet and outlet regions by linear extrapolation (G1 continuity)
        The hub surface is extended one-fourth of the meridional channel arc length at midspan

        """

        # Sorting trick to retrieve order after concatenation
        my_order1 = np.argsort(u)
        my_order2 = np.argsort(my_order1)

        # Define the parameter for each arc (split the vector u in 3 pieces)
        # This is reordering, and we cannot reorder during matching. We retrieve the original order later
        u_inlet = np.sort((u[(u >= 0.00) & (u < 0.25)] - 0.00) / (0.25 - 0.00))
        u_main  = np.sort((u[(u >= 0.25) & (u < 0.75)] - 0.25) / (0.75 - 0.25))
        u_exit  = np.sort((u[(u >= 0.75) & (u <= 1.00)] - 0.75) / (1.00 - 0.75))

        # Get blade arc length at the mean section
        x_func = lambda uu: self.get_meridional_channel_x(uu, 0.50)
        z_func = lambda uu: self.get_meridional_channel_z(uu, 0.50)
        m_func = lambda uu: np.concatenate((x_func(uu), z_func(uu)), axis=0)
        arc_length = get_arc_length(m_func, 0.0 + 1e-6, 1.0 - 1e-6)

        # Tangent line extending the hub surface to the inlet region
        step = 1e-12
        dxdu = np.imag(self.DVs_functions["x_hub"](0 + step * 1j)[0, 0]) / step
        dzdu = np.imag(self.DVs_functions["z_hub"](0 + step * 1j)[0, 0]) / step
        slope_inlet = np.arctan2(dzdu,dxdu)
        x_inlet = self.DVs_control_points["x_hub"][0] + (u_inlet - 1) * np.cos(slope_inlet) * arc_length / 4
        z_inlet = self.DVs_control_points["z_hub"][0] + (u_inlet - 1) * np.sin(slope_inlet) * arc_length / 4

        # Tangent line extending the hub surface to the outlet region
        dxdu = -np.imag(self.DVs_functions["x_hub"](1 - step * 1j)[0, 0]) / step
        dzdu = -np.imag(self.DVs_functions["z_hub"](1 - step * 1j)[0, 0]) / step
        slope_exit = np.arctan2(dzdu,dxdu)
        x_exit = self.DVs_control_points["x_hub"][-1] + u_exit * np.cos(slope_exit) * arc_length / 4
        z_exit = self.DVs_control_points["z_hub"][-1] + u_exit * np.sin(slope_exit) * arc_length / 4

        # Region of the hub surface occupied by the blades
        x_main = self.DVs_functions["x_hub"](u_main).flatten()
        z_main = self.DVs_functions["z_hub"](u_main).flatten()

        # Concatenate the arcs to obtain the extended hub surface (and retrieve the original order)
        x_hub = np.concatenate((x_inlet, x_main, x_exit), axis=0)[my_order2]
        z_hub = np.concatenate((z_inlet, z_main, z_exit), axis=0)[my_order2]
        hub_coordinates = np.asarray((x_hub, z_hub))

        return hub_coordinates


    # ---------------------------------------------------------------------------------------------------------------- #
    # Compute the coordinates of the shroud
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_shroud_surface(self):
        self.u_shroud = np.linspace(0, 1, 100)
        self.shroud_coordinates = self.get_extended_shroud_coordinates(self.u_shroud)

    def get_extended_shroud_coordinates(self, u):

        """ Compute the coordinates of the hub surface in the (x,z) plane

        The shroud surface is is extended to the inlet and outlet regions by linear extrapolation (G1 continuity)
        The shroud surface is extended one-fourth of the meridional channel arc length at midspan

        """

        # Sorting trick to retrieve order after concatenation
        my_order1 = np.argsort(u)
        my_order2 = np.argsort(my_order1)

        # Define the parameter for each arc (split the vector u in 3 pieces)
        # This is reordering, and we cannot reorder during matching. We retrieve the original order later
        u_inlet = np.sort((u[(u >= 0.00) & (u < 0.25)] - 0.00) / (0.25 - 0.00))
        u_main = np.sort((u[(u >= 0.25) & (u < 0.75)] - 0.25) / (0.75 - 0.25))
        u_exit = np.sort((u[(u >= 0.75) & (u <= 1.00)] - 0.75) / (1.00 - 0.75))

        # Get blade arc length at the mean section
        x_func = lambda uu: self.get_meridional_channel_x(uu, 0.50)
        z_func = lambda uu: self.get_meridional_channel_z(uu, 0.50)
        m_func = lambda uu: np.concatenate((x_func(uu), z_func(uu)), axis=0)
        arc_length = get_arc_length(m_func, 0.0 + 1e-6, 1.0 - 1e-6)

        # Tangent line extending the hub surface to the inlet region
        step = 1e-12
        dxdu = np.imag(self.DVs_functions["x_shroud"](0 + step * 1j)[0, 0]) / step
        dzdu = np.imag(self.DVs_functions["z_shroud"](0 + step * 1j)[0, 0]) / step
        slope_inlet = np.arctan2(dzdu, dxdu)
        x_inlet = self.DVs_control_points["x_shroud"][0] + (u_inlet - 1) * np.cos(slope_inlet) * arc_length / 4
        z_inlet = self.DVs_control_points["z_shroud"][0] + (u_inlet - 1) * np.sin(slope_inlet) * arc_length / 4

        # Tangent line extending the hub surface to the outlet region
        dxdu = -np.imag(self.DVs_functions["x_shroud"](1 - step * 1j)[0, 0]) / step
        dzdu = -np.imag(self.DVs_functions["z_shroud"](1 - step * 1j)[0, 0]) / step
        slope_exit = np.arctan2(dzdu, dxdu)
        x_exit = self.DVs_control_points["x_shroud"][-1] + u_exit * np.cos(slope_exit) * arc_length / 4
        z_exit = self.DVs_control_points["z_shroud"][-1] + u_exit * np.sin(slope_exit) * arc_length / 4

        # Region of the hub surface occupied by the blades
        x_main = self.DVs_functions["x_shroud"](u_main).flatten()
        z_main = self.DVs_functions["z_shroud"](u_main).flatten()

        # Concatenate the arcs to obtain the extended hub surface (and retrieve the original order)
        x_hub = np.concatenate((x_inlet, x_main, x_exit), axis=0)[my_order2]
        z_hub = np.concatenate((z_inlet, z_main, z_exit), axis=0)[my_order2]
        shroud_coordinates = np.asarray((x_hub, z_hub))

        return shroud_coordinates
    

    def get_section_thickness_properties(self,section_upper_side,section_lower_side):
    
        """ Compute the chordwise thickness distribution between the upper and lower side curves
        For optimisation purposes, this will allow to impose constraints such as:
        - A limit on or a specific radial distribution of maximum thickness is desired
        - Except for very particular cases, a well designed airfoil only has one global
         maximum, i.e. dTh/dx is monotonic
   
        Parameters
        ----------
        section_upper_side: Array of upper side x-y coordinates
        
        section_lower_side: Array of lower side x-y coordinates
       
        Returns
        -------
        section_thickness_dist : Array of shape (2,points). For each section, its axial thickness distribution
            
        maxt_th                : The maximum thickness radial distribution
       
        lack_of_monotonicity   : A measure of the lack of monotonicity of the axial derivative of the thickness 
                                If lower than or equal to 0, it is monotone
                                Given as a radial distribution


        Author: Ricardo Puente, 09/2020
                r.puente@imperial.ac.uk
        """
        
        # Compute the camber line as the bisector between the lower and the upper side
        camber, upper_tangent,lower_tangent = get_bisectors(section_upper_side,section_lower_side)   
        
        # # DBG plots
        #     # Airfoil surfaces and camber line            
        # plot3fun(section_upper_side.real,'SS',section_lower_side.real,'PS', camber.real,'Camber')
        # ####
        
        # Compute the thickness distribution as the distance between the tangent points
        m = len(camber[0])
        section_thickness_dist = np.zeros((2,m),dtype=complex)
        for i in range(m):
            up = upper_tangent[i]
            lp = lower_tangent[i]
            dv = up-lp
            d = (dv[0]*dv[0]+dv[1]*dv[1])**(0.5)
            section_thickness_dist[1][i] = d
            section_thickness_dist[0][i] = camber[0][i]
             
        # Get maximum
        max_th = max(section_thickness_dist[1])
        
        # Compute derivative of thickness distribution
        thickness_der = get_curve_derivative(camber[0],section_thickness_dist[1])
        
        # Get monotonicity measure
        lack_of_monotonicity = get_monotonicity_measure(camber[0],thickness_der)
        
        # # DBG plots
        #     # Thickness and thickness derivative
        # plotfun_xy(camber[0].real,section_thickness_dist.real,'Thickness')
        # plotfun_xy(camber[0].real,thickness_der.real,'Thickness Derivative')        
        # #############
        
        return section_thickness_dist,max_th,lack_of_monotonicity
    
    def get_blade_thickness_properties(self):
    
        """ Call get_section_thickness_properties for every section to give the radial distribution
            of thickness properties for each control section   
       
        Returns
        -------
        section_thickness_dist    : Axial thickness distribution per section
        
        max_th_dist               : The maximum thickness per section
       
        lack_of_monotonicity_dist : Monotonicity measure per section


        Author: Ricardo Puente, 09/2020
                r.puente@imperial.ac.uk
        """
        
        section_thickness_dist    = []
        max_th_dist               = np.zeros((self.N_SECTIONS,),dtype=complex)
        lack_of_monotonicity_dist = np.zeros((self.N_SECTIONS,),dtype=complex)
        
        r_coor = self.make_radial_distribution()

        for i in range(0,self.N_SECTIONS):
            upper_side,lower_side = split_curve(self.section_coordinates[i])
            section_thickness,max_th_dist[i],lack_of_monotonicity_dist[i] = self.get_section_thickness_properties(upper_side,lower_side)
            section_thickness_dist.append(section_thickness)            
        
        return section_thickness_dist,max_th_dist,lack_of_monotonicity_dist,r_coor
