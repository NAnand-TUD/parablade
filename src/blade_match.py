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

#---------------------------------------------------------------------------------------------#
# Importing general packages
#---------------------------------------------------------------------------------------------#
import sys                                              # Throw exception errors
import os                                               # Get current working directory
import time                                             # Quick time measurement
import pdb                                              # Python debugging tool
import copy                                             # Make a deep copy of an object
from scipy.optimize import *                            # Optimization library
import numpy as np                                      # Scientific computing library
try:
    import matplotlib as mpl                                # Plotting library
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.animation as animation
except:
    pass
from scipy import stats
import warnings

#---------------------------------------------------------------------------------------------#
# Importing user-defined packages
#---------------------------------------------------------------------------------------------#
from blade_3D import Blade3D
from common import printProgress
from config import ReadUserInput, WriteBladeConfigFile


#---------------------------------------------------------------------------------------------#
# Define the blade matching class
#---------------------------------------------------------------------------------------------#
class BladeMatch:


    """ Create a BladeMatch object to match the geometry of an prescribed blade

    BladeMatch can be used to find the set of design variables of a Blade3D object that replicate the geometry of an
    existing turbomachinery blade. This is an important step in shape optimization because it allows to find a
    parametric representation of the baseline geometry. This geometry can then be modified in a systematic way by an
    optimization algorithm to maximize or minimize a suitable performance indicator

    Parameters
    ----------
    IN_file : string
        Path to the configuration file that contains the initial guess for the matching
        The configuration file must specify the field PRESCRIBED_BLADE_FILENAME with the name of the file that contains
        the geometry of the prescribed blade
        The coordinates of the prescribed blade should be specified as tab-separated (index, x, y, z) columns for 3D
        problems or as tab-separated (index, x, y) columns for 2D problems

    coarseness : int
        Factor used to sample one of every ´coarseness´-th point of the prescribed geometry
        This may be useful to reduce the computational time of the blade matching during preliminary tests


    References
    ----------
    TODO add publication information here

    """

    def __init__(self, IN, coarseness=1, plot_options=None):

        # Declare input variables as instance variables
        self.IN                         = IN
        self.NDIM                       = int(self.IN["NDIM"][0])
        self.N_SECTIONS                 = self.IN["N_SECTIONS"][0]
        self.PRESCRIBED_BLADE_FILENAME  = self.IN["PRESCRIBED_BLADE_FILENAME"]
        self.plot_options               = plot_options

        # Create output directory
        os.system("rm -rf output_matching")
        os.mkdir("output_matching")

        # Load prescribed blade coordinates
        if self.NDIM == 2:
            self.coordinates_prescribed = np.loadtxt(self.PRESCRIBED_BLADE_FILENAME, delimiter='\t').transpose()
            self.coordinates_prescribed = self.coordinates_prescribed[:, ::coarseness]
            self.N_points = np.shape(self.coordinates_prescribed)[1]
            self.coordinates_prescribed = np.concatenate((self.coordinates_prescribed, np.zeros((1, self.N_points))))
        elif self.NDIM == 3:
            self.coordinates_prescribed = np.loadtxt(self.PRESCRIBED_BLADE_FILENAME, delimiter='\t').transpose()
            # self.coordinates_prescribed = self.coordinates_prescribed[:, ::coarseness]
            self.coordinates_prescribed = self.coordinates_prescribed[[0, 3, 2, 1], ::coarseness] # Fix SU2 convention
            # R_condition = (self.coordinates_prescribed[[3],:]**2 + self.coordinates_prescribed[[1],:]**2)**(1/2)
            # self.coordinates_prescribed = np.where(R_condition < 0.046, self.coordinates_prescribed, 0*self.coordinates_prescribed)
            self.N_points = np.shape(self.coordinates_prescribed)[1]
        else: raise Exception('The number of dimensions must be 2 or 3')

        # Initialize optimization problem
        self.function_calls = 0
        self.iteration = 0

        # Initialize mismatch indicators
        self.mean_deviation = 0
        self.max_deviation = 0

        # Initialize the matched blade object
        self.blade_matched = Blade3D(self.IN)
        self.blade_matched.initialize_uv_values(Nu=200, Nv=int(self.N_SECTIONS))
        self.blade_matched.make_blade()
        self.u = self.blade_matched.u
        self.v = self.blade_matched.v
        self.coordinates_matched = self.blade_matched.get_surface_coordinates(self.u, self.v)
        self.error_distribution = None


    # ---------------------------------------------------------------------------------------------------------------- #
    # Blade matching main function
    # ---------------------------------------------------------------------------------------------------------------- #
    def match_blade(self, matching_mode='auto'):

        """ Match a blade parametrization to a prescribed blade

        This function contains three matching modes:

        1) Manual mode: matching_mode='auto'

            - Plots the prescribed and fitted blades in an interactive way (read instruction on screen)
            - The geometry of the matched blade can be updated in real time editing the .cfg file
            - The output from this matching mode is used as initial guess or the automatic matching

        2) (u,v) mode: matching_mode='uv'

              Solves N independent optimization problems using the parameters (u,v) as independent variables and finds
              the (u,v) pair that minimizes the error between the prescribed and parametrized blades.
              Each optimization problem is solved using different starting points to find the global minima of the
              error and avoid getting stuck in local minima
              This mode is useful to get (u,v) parametrization of the blade surface points of a mesh
              (both structured and unstructured meshed can be matched)

        3) Automatic mode: matching_mode='DVs'

              Solves a single optimization problem using the design variables as independent variables and finds the
              design variables that minimize the error between the prescribed and parametrized blades
              When this mode is executed the (u,v) parametrization is matched once and then the design variable matching
              starts. The (u,v) matching is re-executed each N iterations of the design variable matching to avoid
              getting stuck in (u,v)-local minima (see callback func)

        """

        if matching_mode == 'manual':

            # Plot the prescribed and matched blades in interactive mode
            self.plot_blade_matching()
            self.do_interactive_matching()

        elif matching_mode == 'uv':

            # Run the (u,v) matching
            self.match_blade_uv()

            # Plot the error lines and save the final solution
            self.plot_blade_matching()
            self.plot_error_distribution()
            self.plot_error_lines()
            self.save_plots()
            self.print_coordinates()
            self.print_config_file()
            plt.show()

        elif matching_mode == 'DVs':

            # Run the (u,v) matching one time first
            self.match_blade_uv()

            # Plot the initial guess in a new figure
            self.plot_blade_matching()
            self.plot_error_distribution()
            self.update_plots()

            # Run the design variable matching | The (u,v) matching is run each N iterations (see callback function)
            self.match_blade_DVs()

            # Plot the error lines and save the final solution
            self.update_plots()
            self.plot_error_lines()
            self.save_plots()
            self.print_coordinates()
            self.print_config_file()
            plt.show()

        else:

            raise Exception("Choose a valid option for matching_mode: 'manual', 'uv', or 'DVs'")


    # ---------------------------------------------------------------------------------------------------------------- #
    # Match the (u,v) parametrization
    # ---------------------------------------------------------------------------------------------------------------- #
    def match_blade_uv(self):

        """ Solve the blade matching problem for each prescibed point using (u,v) as independent variables """

        # Initialize the arrays for the (u,v) parametrization
        print('\n', 'Starting (u,v) parametrization matching...')
        my_u = []
        my_v = []
        h = 1e-5

        # Solve N independent optimization problems
        for i in range(self.N_points):

            # Display the matching progress
            printProgress(i, self.N_points)

            # Start the (u,v) matching from different initial values
            if self.NDIM == 2:
                my_u0 = [0.100, 0.250, 0.500, 0.750, 0.900]
                my_v0 = [0.500, 0.500, 0.500, 0.500, 0.500]
                my_bounds = [[(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)]]
            else:
                my_u0 = [0.100, 0.250, 0.500, 0.750, 0.900, 0.100, 0.250, 0.500, 0.750, 0.900]
                my_v0 = [0.250, 0.250, 0.250, 0.250, 0.250, 0.750, 0.750, 0.750, 0.750, 0.750]
                my_bounds = [[(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
                             [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)]]

            # Add the converged solution from the previous iteration as initial guess
            if self.iteration > 0:
                my_u0.append(self.u[i])
                my_v0.append(self.v[i])
                my_bounds.append([(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)])

            # Initialize lists to save the solution for each starting point
            my_fun = []
            my_x = []

            # Solve the optimization problem
            for k in range(len(my_u0)):

                # Optimization algorithm options
                my_options = {'disp': False,
                              'ftol': 1e-8,
                              #'gtol': 1e-9,
                              #'eps': np.
                              # finfo(np.float64).eps ** (1 / 2),
                              'maxiter': 1000}

                # Solve the optimization problem
                solution = minimize(fun=self.my_objective_function,
                                    x0=np.asarray([my_u0[k], my_v0[k]]),
                                    args=('uv_parametrization', i),
                                    method='SLSQP',   # 'SLSQP' proved to be more robust and faster than 'L-BFGS-B'
                                    jac=None,
                                    # hess=None,
                                    # hessp=None,
                                    bounds=my_bounds[k],
                                    # constraints=None,
                                    callback=None,
                                    options=my_options)

                # Store the current solution
                my_x.append(solution.x)
                my_fun.append(solution.fun)

            # Pick the global minimum
            index = int(np.argmin(my_fun))
            my_u.append(my_x[index][0])
            my_v.append(my_x[index][1])

        # Update the class (u,v) parametrization
        self.u = np.asarray(my_u)
        self.v = np.asarray(my_v)
        self.coordinates_matched = self.blade_matched.get_surface_coordinates(self.u, self.v)
        print('\n', '(u,v) parametrization matching is finished...')


    # ---------------------------------------------------------------------------------------------------------------- #
    # Match design variables
    # ---------------------------------------------------------------------------------------------------------------- #
    def match_blade_DVs(self):

        """ Solve the blade matching problem using the design variables as independent variables """

        # Initialize design variable names and bounds
        print('\n', 'Starting design variables matching...')
        self.DVs_names = self.blade_matched.DVs_names
        # self.DVs_names = self.blade_matched.DVs_names_2D
        # self.DVs_names = self.blade_matched.DVs_names_meridional

        variable_values = []
        for k in self.DVs_names:
            variable_values.extend(self.IN[k])

        # Initialize the design variable initial guess (starting point)
        my_x0 = np.asarray(variable_values)

        # Initialize the design variable bounds
        self.initialize_DVs_bounds()
        my_bounds = self.DVs_bounds

        # Start the iteration and function call counter
        self.function_calls = 0
        self.iteration = 0

        # Optimization algorithm options
        my_options = {'disp': False,
                      'ftol': 1e-8,
                      # 'gtol': 1e-9,
                      # 'eps': np.finfo(np.float64).eps ** (1 / 2),
                      'maxiter': 1000}

        # Solve the optimization problem
        self.solution = minimize(fun=self.my_objective_function,
                                 x0=my_x0,
                                 args='design_variables',
                                 method='SLSQP',
                                 jac=None,
                                 # hess=None,
                                 # hessp=None,
                                 bounds=my_bounds,
                                 # constraints=None,
                                 callback=self.callback_function,
                                 options=my_options)

        self.coordinates_matched = self.blade_matched.get_surface_coordinates(self.u, self.v)
        print('\n', 'Design variable matching is finished...')


    # ---------------------------------------------------------------------------------------------------------------- #
    # Initialize the design variable bounds
    # ---------------------------------------------------------------------------------------------------------------- #
    def initialize_DVs_bounds(self):

        """ Set the bounds for the design variables

            The blade matching problem could be formulated as an unconstrained optimization problem

            However, experience showed that it is necessary to specify some reasonable bounds so that the internal
            iterations do not lead to strange unfeasible geometries that would prevent the optimization to converge

            The optimization bounds set by this function constrain the internal iterations but should be inactive
            constraints when the optimization converges

        """

        # Use the meanline arc-length as a reference length to set the bounds
        self.meanline_length = np.real(self.blade_matched.meanline_length)

        # Initialize an empty list to store the tuples containing the design variable bounds
        self.DVs_bounds = []

        # Iterate over all design variable names and number of control points
        for k in self.DVs_names:
            for i in range(len(self.IN[k])):

                if k in ['x_leading', 'y_leading', 'z_leading', 'x_trailing', 'z_trailing',
                         'x_hub', 'z_hub', 'x_shroud', 'z_shroud']:
                    self.DVs_bounds.append((self.IN[k][i] - self.meanline_length, self.IN[k][i] + self.meanline_length))

                elif k in ['theta_in', 'theta_out', 'stagger']:
                    self.DVs_bounds.append((-89.9, 89.9))

                elif k in ['wedge_in', 'wedge_out']:
                    self.DVs_bounds.append((0.1, 89.9))

                elif k in ['radius_in', 'radius_out']:
                    self.DVs_bounds.append((1e-4, 1.000))

                elif k in ['dist_in', 'dist_out', 'dist_1', 'dist_2', 'dist_3', 'dist_4']:
                    self.DVs_bounds.append((1e-4, 2.000))

                elif k in ['thickness_upper_1', 'thickness_upper_2', 'thickness_upper_3',
                           'thickness_upper_4', 'thickness_upper_5', 'thickness_upper_6',
                           'thickness_lower_1', 'thickness_lower_2', 'thickness_lower_3',
                           'thickness_lower_4', 'thickness_lower_5', 'thickness_lower_6']:
                    self.DVs_bounds.append((1e-4, 1.000))


    # ---------------------------------------------------------------------------------------------------------------- #
    # Define the objective function of the problem
    # ---------------------------------------------------------------------------------------------------------------- #
    def my_objective_function(self, x, optimization_type, i=0):

        """ Evaluate the objective function of the optimization problem

            The objective function is the two-norm of the mismatch between the prescribed and the matched blades

            When optimization_type='uv' the vector of degrees of freedom contains the current (u,v) parametrization.
            In this case the error between prescribed and fitted blades is evaluated at one point at a time

            When optimization_type='DVs' the vector of degrees of freedom contains the design variables
            In this case the error between prescribed and fitted blades is evaluated at every point at once (two-norm)

        """

        if optimization_type == 'uv_parametrization':

            # Use the parameters (u,v) as independent variables
            u = np.asarray([x[0]])
            v = np.asarray([x[1]])

        elif optimization_type == 'design_variables':

            # Use the design variables as independent variables
            variable_values = x
            k_start=0
            for key in self.DVs_names:
                len_value = len(self.IN[key])
                k_end = k_start+len_value
                self.IN[key] = variable_values[k_start:k_end].tolist()
                k_start = k_end

            # Update blade design variables
            self.blade_matched.update_DVs_control_points(self.IN)
            self.blade_matched.make_blade()

            # Get (u,v) parametrization
            u = self.u
            v = self.v

        else:
            raise Exception("Choose a valid option for optimization_type: 'parameter_distribution' or 'design_variables'")


        # Get the coordinates of the matched blade
        coordinates_matched = self.blade_matched.get_surface_coordinates(u, v)

        # Get the coordinates of the prescribed blade
        if optimization_type == 'uv_parametrization':
            coordinates_prescribed = self.coordinates_prescribed[1:, [i]]
        else:
            coordinates_prescribed = self.coordinates_prescribed[1:, :]

        # Ignore the z-coordinate for 2D problems
        if self.NDIM == 2:
            coordinates_matched = coordinates_matched[[0, 1], :]
            coordinates_prescribed = coordinates_prescribed[[0, 1], :]

        # Compute the two-norm of the deviation between prescribed and matched blade
        two_norm_error = np.real(np.sum((coordinates_matched - coordinates_prescribed) ** 2) ** (1 / 2))

        # Compute additional mismatch indicators
        if optimization_type == 'design_variables':
            self.error_distribution = np.real(np.sum((coordinates_matched - coordinates_prescribed) ** 2, axis=0) ** (1 / 2))
            self.mean_deviation = self.error_distribution.mean()
            self.mean_deviation_rel = self.mean_deviation / self.meanline_length * 100
            self.max_deviation = np.real(np.max(np.sum((coordinates_matched-coordinates_prescribed)**2,axis=0)**(1/2)))
            self.max_deviation_rel = self.max_deviation / self.meanline_length * 100

        # Update number of function calls
        self.function_calls = self.function_calls + 1

        return two_norm_error


    # ---------------------------------------------------------------------------------------------------------------- #
    # Optimization callback function
    # ---------------------------------------------------------------------------------------------------------------- #
    def callback_function(self, x):

        """ Display optimization progress and update the plots

            In addition callback_function() performs the (u,v) parametrization matching every N-th iteration

        """

        # Update iteration number
        self.iteration = self.iteration + 1

        # Print header each N iterations
        N = 10
        if (self.iteration == 1) or ((self.iteration-1) % N == 0):
            print('\n Iteration \t Mean deviation (m) \t Maximum deviation (m) \t Mean deviation (%) \t Maximum deviation (%)' )

        # Print the optimization status
        print('{:>10} \t {:>18.6f} \t {:>21.6f} \t {:>18.6f} \t {:>21.6f}'.format(
            self.iteration, self.mean_deviation, self.max_deviation, self.mean_deviation_rel, self.max_deviation_rel))

        # Update the plot each iteration
        if self.iteration % 1 == 0:
            self.update_plots()

        # Update the (u,v)-parametrization each N iterations
        if self.iteration % N == 0:
            self.match_blade_uv()

        # Update the matched coordinates
        self.coordinates_matched = self.blade_matched.get_surface_coordinates(self.u, self.v)

        # Print the matching status and optimization progress
        self.print_coordinates()
        self.print_config_file()
        self.print_optimization_progress()


    # ---------------------------------------------------------------------------------------------------------------- #
    # Save the current blade matching
    # ---------------------------------------------------------------------------------------------------------------- #
    def print_coordinates(self, filename='matched_coordinates', path=os.getcwd()):

        """ Print the coordinates of the prescribed and matched blades in a .csv file """

        # Print the matched blade coordinates
        full_path = path + '/output_matching/'
        file = open(full_path + filename + '.csv', 'w')

        if self.NDIM == 2:
            file.write("\"index\",\t\"x_prescribed\",\t\"y_prescribed\",\t\"x_match\",\t\"y_match\",\t\"u\",\t\"v\"\n")
            for i in range(self.N_points):
                file.write('%i,,\t%+.15e,,\t%+.15e,,\t%+.15e,,\t%+.15e,,\t%.15f,,\t%.15f\n' %
                           (self.coordinates_prescribed[0, i],          # Mesh point index
                            self.coordinates_prescribed[1, i],          # Prescribed x coordinate (axial)
                            self.coordinates_prescribed[2, i],          # Prescribed y coordinate (tangential)
                            np.real(self.coordinates_matched[0, i]),    # Matched x-coordinate (axial)
                            np.real(self.coordinates_matched[1, i]),    # Matched y-coordinate (tangential)
                            self.u[i], self.v[i]))                      # (u,v) parametrization

        elif self.NDIM == 3:
            file.write("\"index\",\t\"x_prescribed\",\t\"y_prescribed\",\t\"z_prescribed\",\t\"x_match\",\t\"y_match\",\t\"z_match\",\t\"u\",\t\"v\"\n")
            for i in range(self.N_points):
                file.write('%i,\t%+.15e,\t%+.15e,\t%+.15e,\t%+.15e,\t%+.15e,\t%+.15e,\t%.15f,\t%.15f\n' %
                           (self.coordinates_prescribed[0, i],          # Mesh point index
                            self.coordinates_prescribed[3, i],          # SU2 x coordinate (radial)
                            self.coordinates_prescribed[2, i],          # SU2 y coordinate (tangential)
                            self.coordinates_prescribed[1, i],          # SU2 z coordinate (axial)
                            np.real(self.coordinates_matched[2, i]),    # SU2 x coordinate (radial)
                            np.real(self.coordinates_matched[1, i]),    # SU2 y coordinate (tangential)
                            np.real(self.coordinates_matched[0, i]),    # SU2 z coordinate (axial)
                            self.u[i], self.v[i]))                      # SU2 z coordinate (axial)

        else:
            raise Exception('The number of dimensions must be 2 or 3')

        file.close()


    def print_config_file(self, filename="matched_parametrization", path=os.getcwd()):

        """ Print a configuration .cfg file for the current set of design variables """

        full_path = path + '/output_matching/'
        file = open(full_path + filename + '.cfg', 'w')
        WriteBladeConfigFile(file, self.IN)
        file.close()


    def print_optimization_progress(self, filename="optimization_progress", path=os.getcwd()):

        """ Print the convergence history of the design variable optimization """

        # Open the file in "append" mode
        full_path = path + '/output_matching/'
        file = open(full_path + filename + '.txt', 'a')

        # Print the header for the first iteration
        if self.iteration == 1:
            header = 'Iteration \t Mean deviation (m) \t Maximum deviation (m) \t Mean deviation (%) \t Maximum deviation (%)\n'
            file.write(header)

        # Print optimization progress
        file.write('%i \t %.6f \t %.6f \t %.6f \t %.6f\n' %
                            (self.iteration,
                             self.mean_deviation,
                             self.max_deviation,
                             self.mean_deviation_rel,
                             self.max_deviation_rel))

        file.close()


    # ---------------------------------------------------------------------------------------------------------------- #
    # Make plots for the matching
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_blade_matching(self):

        """ Plot of the prescribed and the matched blades

            The figures plotted can be controlled with the "options" dictionary
            The plotting possibilities include:

                1) 2D view of the (z,y)-plane (airfoil sections)
                2) 2D view of the (z,R)-plane (meridional channel)
                3) 2D view of the (y,x)-plane (front view of the blade)
                4) 3D view of the blade

        """

        # Prescribed coordinates
        x_prescribed = self.coordinates_prescribed[1, :]
        y_prescribed = self.coordinates_prescribed[2, :]
        z_prescribed = self.coordinates_prescribed[3, :]
        R_prescribed = np.sqrt(y_prescribed ** 2 + z_prescribed ** 2)

        # Matched coordinates
        x_matched = np.real(self.coordinates_matched[0, :])
        y_matched = np.real(self.coordinates_matched[1, :])
        z_matched = np.real(self.coordinates_matched[2, :])
        R_matched = np.sqrt(y_matched ** 2 + z_matched ** 2)

        # Get minimum axes values
        x_max = np.max(x_prescribed)
        y_max = np.max(y_prescribed)
        z_max = np.max(z_prescribed)
        R_max = np.max(R_prescribed)

        # Get maximum axes values
        x_min = np.min(x_prescribed)
        y_min = np.min(y_prescribed)
        z_min = np.min(z_prescribed)
        R_min = np.min(R_prescribed)

        # Get the arithmetic centre of the blade
        x_mid = (x_min + x_max) / 2
        y_mid = (y_min + y_max) / 2
        z_mid = (z_min + z_max) / 2
        R_mid = (R_min + R_max) / 2

        # Activate interactive mode plotting for the animations (not necessary)
        # plt.ion()

        if self.plot_options['view_xy'] == 'yes':

            # Plot the prescribed and the matched blades
            self.figure_1, self.axes_1 = plt.subplots()
            self.axes_1.set_aspect(1.00)
            self.axes_1.set_xlabel('$x$ - axis', fontsize=12, color='k', labelpad=12)
            self.axes_1.set_ylabel('$y$ - axis', fontsize=12, color='k', labelpad=12)
            # self.axes_1.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            # self.axes_1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            for t in self.axes_1.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in self.axes_1.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            plt.tight_layout(pad=1.0, w_pad=None, h_pad=None)

            # Set the axes for the plot
            L = np.max((x_max - x_min, y_max - y_min)) / 2
            self.axes_1.set_xlim([x_mid - 1.25 * L, x_mid + 1.25 * L])
            self.axes_1.set_ylim([y_mid - 1.25 * L, y_mid + 1.25 * L])

            # Plot prescribed coordinates
            self.points_1p, = self.axes_1.plot(x_prescribed, y_prescribed)
            self.points_1p.set_marker("o")
            self.points_1p.set_markersize(2.5)
            self.points_1p.set_markeredgewidth(0.5)
            self.points_1p.set_markeredgecolor("k")
            self.points_1p.set_markerfacecolor("w")
            self.points_1p.set_linestyle(" ")
            self.points_1p.set_color("k")
            self.points_1p.set_linewidth(0.50)
            self.points_1p.set_label('Blade prescribed')

            # Plot matched coordinates
            self.points_1m, = self.axes_1.plot(x_matched, y_matched)
            self.points_1m.set_marker("o")
            self.points_1m.set_markersize(2.5)
            self.points_1m.set_markeredgewidth(0.5)
            self.points_1m.set_markeredgecolor("b")
            self.points_1m.set_markerfacecolor("w")
            self.points_1m.set_linestyle(" ")
            self.points_1m.set_color("k")
            self.points_1m.set_linewidth(0.50)
            self.points_1m.set_label('Blade matched')
            

        if self.plot_options['view_xR'] == 'yes':

            # Plot the prescribed and the matched blades
            self.figure_2, self.axes_2 = plt.subplots()
            self.axes_2.set_aspect(1.00)
            self.axes_2.set_xlabel('$x$ - axis', fontsize=12, color='k', labelpad=12)
            self.axes_2.set_ylabel('$R$ - axis', fontsize=12, color='k', labelpad=12)
            # self.axes_2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            # self.axes_2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            for t in self.axes_2.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in self.axes_2.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            plt.tight_layout(pad=1.0, w_pad=None, h_pad=None)

            # Set the axes for the plot
            L = np.max((x_max - x_min, R_max - R_min)) / 2
            self.axes_2.set_xlim([x_mid - 1.25 * L, x_mid + 1.25 * L])
            self.axes_2.set_ylim([R_mid - 1.25 * L, R_mid + 1.25 * L])

            # Plot prescribed coordinates
            self.points_2p, = self.axes_2.plot(x_prescribed, R_prescribed)
            self.points_2p.set_marker("o")
            self.points_2p.set_markersize(3.0)
            self.points_2p.set_markeredgewidth(0.5)
            self.points_2p.set_markeredgecolor("k")
            self.points_2p.set_markerfacecolor("w")
            self.points_2p.set_linestyle(" ")
            self.points_2p.set_color("k")
            self.points_2p.set_linewidth(0.50)
            self.points_2p.set_label('Blade prescribed')

            # Plot matched coordinates
            self.points_2m, = self.axes_2.plot(x_matched, R_matched)
            self.points_2m.set_marker("o")
            self.points_2m.set_markersize(3.0)
            self.points_2m.set_markeredgewidth(0.5)
            self.points_2m.set_markeredgecolor("b")
            self.points_2m.set_markerfacecolor("w")
            self.points_2m.set_linestyle(" ")
            self.points_2m.set_color("k")
            self.points_2m.set_linewidth(0.50)
            self.points_2m.set_label('Blade matched')


        if self.plot_options['view_yz'] == 'yes':

            # Plot the prescribed and the matched blades
            self.figure_3, self.axes_3 = plt.subplots()
            self.axes_3.set_aspect(1.00)
            self.axes_3.set_xlabel('$y$ - axis', fontsize=12, color='k', labelpad=12)
            self.axes_3.set_ylabel('$z$ - axis', fontsize=12, color='k', labelpad=12)
            # self.axes_3.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            # self.axes_3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            for t in self.axes_3.xaxis.get_major_ticks(): t.label.set_fontsize(12)
            for t in self.axes_3.yaxis.get_major_ticks(): t.label.set_fontsize(12)
            plt.tight_layout(pad=1.0, w_pad=None, h_pad=None)

            # Set the axes for the plot
            L = np.max((y_max - y_min, z_max - z_min)) / 2
            self.axes_3.set_xlim([y_mid - 1.25 * L, y_mid + 1.25 * L])
            self.axes_3.set_ylim([z_mid - 1.25 * L, z_mid + 1.25 * L])

            # Plot prescribed coordinates
            self.points_3p, = self.axes_3.plot(y_prescribed, z_prescribed)
            self.points_3p.set_marker("o")
            self.points_3p.set_markersize(3.0)
            self.points_3p.set_markeredgewidth(0.5)
            self.points_3p.set_markeredgecolor("k")
            self.points_3p.set_markerfacecolor("w")
            self.points_3p.set_linestyle(" ")
            self.points_3p.set_color("k")
            self.points_3p.set_linewidth(0.50)
            self.points_3p.set_label('Blade prescribed')

            # Plot matched coordinates
            self.points_3m, = self.axes_3.plot(y_matched, z_matched)
            self.points_3m.set_marker("o")
            self.points_3m.set_markersize(3.0)
            self.points_3m.set_markeredgewidth(0.5)
            self.points_3m.set_markeredgecolor("b")
            self.points_3m.set_markerfacecolor("w")
            self.points_3m.set_linestyle(" ")
            self.points_3m.set_color("k")
            self.points_3m.set_linewidth(0.50)
            self.points_3m.set_label('Blade matched')


        if self.plot_options['view_3D'] == 'yes':

            # Prepare the plot
            self.figure_4 = plt.figure()
            self.axes_4 = Axes3D(self.figure_4)
            self.axes_4.view_init(azim=-55, elev=25)
            self.axes_4.grid(False)
            self.axes_4.xaxis.pane.set_edgecolor('black')
            self.axes_4.yaxis.pane.set_edgecolor('black')
            self.axes_4.zaxis.pane.set_edgecolor('black')
            self.axes_4.xaxis.pane.fill = False
            self.axes_4.yaxis.pane.fill = False
            self.axes_4.zaxis.pane.fill = False
            fontsize = 12
            self.axes_4.set_xlabel('$x$ - axis', fontsize=12, color='k', labelpad=12)
            self.axes_4.set_ylabel('$y$ - axis', fontsize=12, color='k', labelpad=12)
            self.axes_4.set_zlabel('$z$ - axis', fontsize=12, color='k', labelpad=12)
            for t in self.axes_4.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
            for t in self.axes_4.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
            for t in self.axes_4.zaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
            # self.axes_4.set_xticks([])
            # self.axes_4.set_yticks([])
            # self.axes_4.set_zticks([])
            # self.axes_4.axis('off')

            # Set the axes for the plot
            L = np.max((x_max - x_min, y_max - y_min, z_max - z_min)) / 2
            self.axes_4.set_xlim3d(x_mid - 1.1 * L, x_mid + 1.1 * L)
            self.axes_4.set_ylim3d(y_mid - 1.1 * L, y_mid + 1.1 * L)
            self.axes_4.set_zlim3d(z_mid - 1.1 * L, z_mid + 1.1 * L)

            # Plot prescribed coordinates
            self.points_4p, = self.axes_4.plot(x_prescribed, y_prescribed, z_prescribed)
            self.points_4p.set_marker("o")
            self.points_4p.set_markersize(3)
            self.points_4p.set_markeredgewidth(0.5)
            self.points_4p.set_markeredgecolor("k")
            self.points_4p.set_markerfacecolor("w")
            self.points_4p.set_linestyle(" ")
            self.points_4p.set_color("k")
            self.points_4p.set_linewidth(0.50)
            self.points_4p.set_label('Blade prescribed')

            # Plot matched coordinates
            self.points_4m, = self.axes_4.plot(x_matched, y_matched, z_matched)
            self.points_4m.set_marker("o")
            self.points_4m.set_markersize(3)
            self.points_4m.set_markeredgewidth(0.5)
            self.points_4m.set_markeredgecolor("b")
            self.points_4m.set_markerfacecolor("w")
            self.points_4m.set_linestyle(" ")
            self.points_4m.set_color("k")
            self.points_4m.set_linewidth(0.50)
            self.points_4m.set_label('Blade matched')


    def plot_error_distribution(self):

        """ Plot the distribution of the fitting error between prescribed and matched blades"""

        if self.plot_options['error_distribution'] == 'yes':

            # Prepare figure
            self.figure_5, self.axes_5 = plt.subplots()
            self.axes_5.set_xlim(0, 1)
            self.axes_5.set_xlabel('Matching error distribution (mm)', fontsize=11, color='k', labelpad=12)
            self.axes_5.set_xscale('linear')
            self.axes_5.set_ylim(0, 1.25)
            self.axes_5.set_yticks([])
            self.axes_5.set_yticklabels([])
            for t in self.axes_5.yaxis.get_major_ticks(): t.label.set_fontsize(10)

            # Make error distribution plot
            # if self.NDIM == 2: error = np.sum((self.coordinates_prescribed[1:-1, :] - self.coordinates_matched[0:-1,:])**2, axis=0)**(1/2)
            # else: error  = np.sum((self.coordinates_prescribed[1:, :] - self.coordinates_matched)**2, axis=0)**(1/2)
            if self.error_distribution is None: self.error_distribution = np.linspace(0,1,200)
            error = np.real(self.error_distribution)*1000
            error_mean = error.mean()
            kernel = stats.gaussian_kde(error)
            x_values = np.linspace(0, np.amax(error), 200)

            # Plot error distribution
            self.error_plot, = self.axes_5.plot(x_values, kernel(x_values)/np.amax(kernel(x_values)))
            self.error_plot.set_marker(" ")
            self.error_plot.set_markersize(3)
            self.error_plot.set_markeredgewidth(0.5)
            self.error_plot.set_markeredgecolor("b")
            self.error_plot.set_markerfacecolor("w")
            self.error_plot.set_linestyle("-")
            self.error_plot.set_color("b")
            self.error_plot.set_linewidth(1.00)
            self.error_plot.set_label('Error distribution')

            # Plot mean error line
            self.error_mean, = self.axes_5.plot([error_mean, error_mean], [0, kernel(error_mean)/np.amax(kernel(x_values))])
            self.error_mean.set_marker(" ")
            self.error_mean.set_markersize(3)
            self.error_mean.set_markeredgewidth(0.5)
            self.error_mean.set_markeredgecolor("b")
            self.error_mean.set_markerfacecolor("w")
            self.error_mean.set_linestyle("-")
            self.error_mean.set_color("k")
            self.error_mean.set_linewidth(1.25)
            self.error_mean.set_label('Mean error')

            # Manufacturing tolerance reference
            manufacturing_line, = self.axes_5.plot([0.05, 0.05], [0, 1.25], 'r', linewidth=1)
            manufacturing_line.set_label('Manufacturing tolerance')

            # Create legend
            self.axes_5.legend(ncol=1, loc='upper right', fontsize=10, edgecolor='k', framealpha=1.0)

            # Adjust PAD
            plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)



    def plot_error_lines(self):

        """ Draw error lines between the prescribed and the matched blades """

        # Prescribed coordinates
        x_prescribed = self.coordinates_prescribed[1, :]
        y_prescribed = self.coordinates_prescribed[2, :]
        z_prescribed = self.coordinates_prescribed[3, :]
        R_prescribed = np.sqrt(y_prescribed ** 2 + z_prescribed ** 2)

        # Matched coordinates
        x_matched = np.real(self.coordinates_matched[0, :])
        y_matched = np.real(self.coordinates_matched[1, :])
        z_matched = np.real(self.coordinates_matched[2, :])
        R_matched = np.sqrt(y_matched ** 2 + z_matched ** 2)

        # Plot the error lines
        if self.plot_options['view_xy'] == 'yes':
            for i in range(self.N_points):
                error_lines_1, = self.axes_1.plot([x_prescribed[i], x_matched[i]], [y_prescribed[i], y_matched[i]])
                error_lines_1.set_color("r")
                error_lines_1.set_linewidth(0.50)

        if self.plot_options['view_xR'] == 'yes':
            for i in range(self.N_points):
                error_lines_2, = self.axes_2.plot([x_prescribed[i], x_matched[i]], [R_prescribed[i], R_matched[i]])
                error_lines_2.set_color("r")
                error_lines_2.set_linewidth(0.50)

        if self.plot_options['view_yz'] == 'yes':
            for i in range(self.N_points):
                error_lines_3, = self.axes_3.plot([y_prescribed[i], y_matched[i]], [z_prescribed[i], z_matched[i]])
                error_lines_3.set_color("r")
                error_lines_3.set_linewidth(0.50)

        if self.plot_options['view_3D'] == 'yes':
            for i in range(self.N_points):
                error_lines_4, = self.axes_4.plot([x_prescribed[i], x_matched[i]], [y_prescribed[i], y_matched[i]], [z_prescribed[i], z_matched[i]])
                error_lines_4.set_color("r")
                error_lines_4.set_linewidth(0.50)


    def update_plots(self):

        """ Update the plot with the current matched geometry """

        # Rename coordinates
        x_matched = np.real(self.coordinates_matched[0, :])
        y_matched = np.real(self.coordinates_matched[1, :])
        z_matched = np.real(self.coordinates_matched[2, :])
        R_matched = np.sqrt(y_matched ** 2 + z_matched ** 2)

        # Update the plot coordinates and re-draw the plot
        if self.plot_options['view_xy'] == 'yes':
            self.points_1m.set_xdata(x_matched)
            self.points_1m.set_ydata(y_matched)
            plt.pause(0.001)                       # Matplotlib needs a moment before re-drawing
            self.figure_1.canvas.draw()
            self.figure_1.canvas.flush_events()

        if self.plot_options['view_xR'] == 'yes':
            self.points_2m.set_xdata(x_matched)
            self.points_2m.set_ydata(R_matched)
            plt.pause(0.001)
            self.figure_2.canvas.draw()
            self.figure_2.canvas.flush_events()

        if self.plot_options['view_yz'] == 'yes':
            self.points_3m.set_xdata(y_matched)
            self.points_3m.set_ydata(z_matched)
            plt.pause(0.001)
            self.figure_3.canvas.draw()
            self.figure_3.canvas.flush_events()

        if self.plot_options['view_3D'] == 'yes':
            self.points_4m.set_xdata(x_matched)
            self.points_4m.set_ydata(y_matched)
            self.points_4m.set_3d_properties(z_matched)
            plt.pause(0.001)
            self.figure_4.canvas.draw()
            self.figure_4.canvas.flush_events()

        if self.plot_options['error_distribution'] == 'yes':
            # if self.NDIM == 2: error = np.sum((self.coordinates_prescribed[1:-1, :] - self.coordinates_matched[0:-1,:])**2, axis=0)**(1/2)
            # else: error  = np.sum((self.coordinates_prescribed[1:, :] - self.coordinates_matched)**2, axis=0)**(1/2)
            error = np.real(self.error_distribution) * 1000
            error_mean = error.mean()
            kernel = stats.gaussian_kde(error)
            x_values = np.linspace(0, np.amax(error), 200)
            self.error_plot.set_xdata(x_values)
            self.error_plot.set_ydata(kernel(x_values)/np.amax(kernel(x_values)))
            self.error_mean.set_xdata([error_mean, error_mean])
            self.error_mean.set_ydata([0, kernel(error_mean)/np.amax(kernel(x_values))])
            plt.pause(0.001)
            self.figure_5.canvas.draw()
            self.figure_5.canvas.flush_events()


    def save_plots(self):

        os.mkdir("output_matching/figures/")

        if self.plot_options['view_xy'] == 'yes':
            self.figure_1.savefig('output_matching/figures/view_xy.pdf', bbox_inches='tight')
        if self.plot_options['view_xR'] == 'yes':
            self.figure_2.savefig('output_matching/figures/view_xR.pdf', bbox_inches='tight')
        if self.plot_options['view_yz'] == 'yes':
            self.figure_3.savefig('output_matching/figures/view_yz.pdf', bbox_inches='tight')
        if self.plot_options['view_3D'] == 'yes':
            self.figure_4.savefig('output_matching/figures/view_3D.pdf', bbox_inches='tight')
        if self.plot_options['error_distribution'] == 'yes':
            self.figure_5.savefig('output_matching/figures/error_distribution.pdf', bbox_inches='tight')


    def do_interactive_matching(self):

        """ Create an iteractive plot that updates the geometry in real time when the .cfg file is changed """

        print("\n")
        print("Opening the interactive plot...\n")
        print("Instructions: \n")
        print("\t - Edit the .cfg file to modify the geometry in real time")
        print("\t - Try to get a good matching manually by trial and error")
        print("\t - The manual match will be used as initial guess for the optimization")
        print("\t - Close the figures when you are done to continue execution")
        print("\n")

        if self.plot_options['view_xy'] == 'yes':
            # This calls the function animate with a refresh time of "interval" miliseconds
            ani_1 = animation.FuncAnimation(self.figure_1, self.animate, interval=500)

        if self.plot_options['view_xR'] == 'yes':
            ani_2 = animation.FuncAnimation(self.figure_2, self.animate, interval=500)

        if self.plot_options['view_yz'] == 'yes':
            ani_3 = animation.FuncAnimation(self.figure_3, self.animate, interval=500)

        if self.plot_options['view_3D'] == 'yes':
            ani_4 = animation.FuncAnimation(self.figure_4, self.animate, interval=500)

        plt.show()
        print("Closing the interactive plot...")


    def animate(self, i):

        """ Update the geometry of the interactive plot """

        # print("Interactive mode is running...")

        # Read the .cfg file and update the geometry
        while True:
            try:
                self.IN = ReadUserInput(self.IN['Config_Path'])
                self.blade_matched = Blade3D(self.IN)
                self.blade_matched.make_blade()
                self.coordinates_matched = self.blade_matched.get_surface_coordinates(self.u, self.v)
                break
            except:
                print("The .cfg file could not be loaded, trying again...")
                time.sleep(1.5)     # Wait a moment in case the .cfg file is not found

        # Rename coordinates
        x_matched = np.real(self.coordinates_matched[0, :])
        y_matched = np.real(self.coordinates_matched[1, :])
        z_matched = np.real(self.coordinates_matched[2, :])
        R_matched = np.sqrt(y_matched ** 2 + z_matched ** 2)

        # Update the plot coordinates
        if self.plot_options['view_xy'] == 'yes':
            self.points_1m.set_xdata(x_matched)
            self.points_1m.set_ydata(y_matched)

        if self.plot_options['view_xR'] == 'yes':
            self.points_2m.set_xdata(x_matched)
            self.points_2m.set_ydata(R_matched)

        if self.plot_options['view_yz'] == 'yes':
            self.points_3m.set_xdata(y_matched)
            self.points_3m.set_ydata(z_matched)

        if self.plot_options['view_3D'] == 'yes':
            self.points_4m.set_xdata(x_matched)
            self.points_4m.set_ydata(y_matched)
            self.points_4m.set_3d_properties(z_matched)

