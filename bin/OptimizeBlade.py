#!/usr/bin/env python3
###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

############################### FILE NAME: OptimizeBlade.py ###################################
#=============================================================================================#
# author: Pablo Garrido, Roberto, Nitish Anand                                                |
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
from scipy.optimize import *
import numpy as np
import os
import pdb
import sys
import time
from copy import deepcopy
import shutil
import errno, subprocess
import pickle

# Imports that might be problematic to run on clusters
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
except:
    pass

#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common import *
from config import *
from blade_3D import Blade3D
from blade_plot import BladePlot
from blade_match import BladeMatch
from blade_output import BladeOutput

#---------------------------------------------------------------------------------------------#
# Print ParaBlade Banner
#---------------------------------------------------------------------------------------------#
PrintBanner()


#---------------------------------------------------------------------------------------------#
## PROCESS: Initializing the DIR and loading the file
#---------------------------------------------------------------------------------------------#
DIR = os.getcwd() + '/'

try:
    INFile = DIR + sys.argv[-1]
except:
    INFile = DIR + 'Optimization.cfg'
try:
    IN = ReadUserInput(INFile)
except:
    raise Exception('\n\n\n''Something went wrong when reading the configuration file,exiting the program...'
                    '\n\nTo call MakeBlade.py from terminal type:'
                    '\n\tMakeBlade.py <configuration file name>')


#---------------------------------------------------------------------------------------------#
## PROCESS: Shape Optimization Class
#---------------------------------------------------------------------------------------------#

class ShapeOptimization:
    # List of design variables to be optimized.
    # It is same as that in blade_3D.py

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

    # Array of variable names for the two-dimensional blade sections (connecting arcs parametrization)
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

    # Global Design Step counter
    DesignStep      = 0

    # Restart design step counter
    DesignStep_Restart = 0

    # MPI run commands
    MPI_COMMAND     = "mpirun -n %i %s %s"

    # Default config file names for Optimization
    ADJ_CONFIG      = "config_CFD_AD.cfg"
    CFD_CONFIG      = "config_CFD.cfg"
    DEF_CONFIG      = "config_DEF.cfg"
    CAD_CONFIG      = "BladeOpti.cfg"

    # Fodler Strucutre for Optimization
    DESIGN_FOLDER       = DIR               + 'DESIGN/DSN_%02i/'
    DEFORM_FOLDER       = DESIGN_FOLDER     + '/DEFORM/'
    DIRECT_FOLDER       = DESIGN_FOLDER     + '/DIRECT/'
    CAD_FOLDER          = DESIGN_FOLDER     + '/CAD/'
    SENS_OBJ_FOLDER     = DESIGN_FOLDER     + '/ADJOINT_%s/'

    # Folder Structure for FD
    FINDIFF_FOLDER      = DIR               + 'FIND_DIFF/'
    DIRECT_FOLDER_FD    = FINDIFF_FOLDER    + '/DIRECT/'
    DEFORM_FOLDER_FD    = FINDIFF_FOLDER    + '/DEFORM/'

    # Slack notifications
    SLACK_NOTIF     = False

    # Optimization History
    OPT_HIST_FILE   = None
    SENS_HIST       = None

    # Pickle files (restoring)
    OPTSET_PCKL     = None                                                  # Contains variable names, values and parameter
    OBJ_PCKL        = None                                                  # Contains objective function value
    CONS_PCKL       = None                                                  # Contains constraint value
    OBJ_SENS_PCKL   = None                                                  # Contains objective function sensitivities
    CONS_SENS_PCKL  = None                                                  # Contains constraint sensitivities

    # Most similar mesh search
    # Searches (in Optimization mode) for the most similar mesh to perform mesh deformation
    MESH_SEARCH     = False

    # Run in DEBUG MODE
    DEBUG           = False

    def __init__(self,IN):
        self.IN             = IN                                            # Optimization Config File name

        self.SU2_CONFIG     = self.IN['SU2_CONFIG_FILE']                    # SU2 Config File name
        self.CAD_CONFIG     = self.IN['BLADE_PAR_FILE']                     # Fitter Blade parameterization file name
        self.MATCH_BLADE    = self.IN['BLADE_MATCH_FILE']                   # Match blade file name

        self.OPERATION_TYPE = self.IN['OPERATION_TYPE']                     # Operation type defined by the user
        self.NDIM           = self.IN['NDIM'][0]                            # Dimentionality of the problem
        self.NZONE          = int(self.IN['NZONE'][0])                      # Number of the zones of the problem
        
        if self.NZONE == 1:
            self.DESIGN_ZONE        = 1
            self.HIST_FILENAME      = 'history'
            self.SENS_FILENAME      = 'surface_sens.dat'
            self.INDEX              = {'ENTROPY_GENERATION': 3,             
                                    'KINETIC_ENERGY_LOSS': 2,
                                    'TOTAL_PRESSURE_LOSS': 1,
                                    'FLOW_ANGLE_OUT': 7,
                                    'DRAG': 2,
                                    'LIFT': 1,
                                    'NORMAL_GRID_VEL': 2, }
        else:
            self.DESIGN_ZONE        = int(self.IN['DESIGN_ZONE'][0])            
            self.HIST_FILENAME      = 'history_' + str(self.DESIGN_ZONE - 1)
            self.SENS_FILENAME      = 'surface_sens_' + str(self.DESIGN_ZONE - 1) + '.dat'
            self.INDEX              = {'ENTROPY_GENERATION': 48,           
                                    'KINETIC_ENERGY_LOSS': 47,
                                    'TOTAL_PRESSURE_LOSS': 46,
                                    'FLOW_ANGLE_OUT': 52, }

        try:
            self.N_CORES = int(self.IN['N_PROCESSORS'][0])                        # User-defined no. of processors
        except:
            self.N_CORES = 2                                                      # Default no. of processors

        # Initialize and perform checks on run options
        self.runtime_options()

        # Initialize design variables
        self.initialize_DVs()

        # Initialise variables
        self.RefOBJ     = None                                              # Reference Obj Function value
        self.RefSens    = None                                              # Reference sensitivities value
        self.Module     = None                                              # Initialize variable used in the optimiser
        self.RefCons    = None

    def start_process(self):
        # This should ideally be the only function to be called from outside.
        # Outputs for user information
        print("--------------------------------------------------")
        print("                      FILES                       ")
        print("SU2 Configuration file   : %s" % self.SU2_CONFIG)
        print("CAD Configuration file   : %s" % self.CAD_CONFIG)
        print("Blade match file         : %s" % self.MATCH_BLADE)
        print("--------------------------------------------------")
        print("                PROBLEM SETUP                     ")
        print("Operation type           : %s" % self.OPERATION_TYPE)
        print("--------------------------------------------------")

        if self.OPERATION_TYPE == "SHAPE_OPTIMIZATION":
            print("Running Shape Optimization...")
            self.store_ref_blade()
            self.run_optimization()

        elif self.OPERATION_TYPE == "DISCRETE_ADJOINT":
            print("Running discrete adjoint...")
            self.store_ref_blade()
            self.run_disc_adj_validation()

        elif self.OPERATION_TYPE == "FINITE_DIFFERENCE":
            print("Running finite-differences...")
            self.SAVE_VOL_GRIDS = True
            self.store_ref_blade()
            if self.SKIP_REF_CFD == 'YES':
                print("*****************************************************************************")
                print("*** Warning! :: FIND_DIFF_ORDER has been forced to 1 as SkipRefCFD = YES  ***")
                print("*****************************************************************************")
                time.sleep(1)
                IN['FIND_DIFF_ORDER'] = [1]
            self.run_finite_diff_validation()

        elif self.OPERATION_TYPE == "PLOT_VALIDATION":
            print('Running Plot Validation...')
            self.SHOW_PLOT       = True
            self.run_plot_validation()

        elif self.OPERATION_TYPE == "PLOT_OPT_HISTORY":
            print('Plot Opt History')

        elif self.OPERATION_TYPE == "PICKLEIZE":
            print("Starting pickleizing ...")
            # Design step to pickleize
            design_step_pickle = 1
            self.store_ref_blade()
            self.pickleize(design_step_pickle)

        else:
            print("Invalid Input :: %s"%(self.OPERATION_TYPE))

    def run_optimization(self):
        """
        Function to run shape optimisation case
        """
        # Check if DESIGN Folder Exists
        if self.RESTART == 'NO':
            try:
                os.system('rm -rf DESIGN')
            except:
                print('No DESIGN folder found')
            os.mkdir('DESIGN')

            # Create Optimization History file
            self.OPT_HIST_FILE = open('./DESIGN/Optimization_history.dat','w')
            self.OPT_HIST_FILE.write("VARIABLES=\"DESIGN_STEP\",\"OBJECTIVE\",\"CONSTRAIN\"\n")

        elif self.RESTART == 'YES':
            self.OPT_HIST_FILE = open("./DESIGN/Optimization_history.dat","a")

        # Initialize design variables from the fitter CAD-parametrization config file.
        self.variable_values    = [None] * self.N_variables         # List that contains initial variable values
        len_values              = [None] * self.N_variables         # List containing length for each DV

        for k in range(self.N_variables):
            self.variable_values[k] = (self.CAD_CONFIG_IN[self.variable_names[k]])
            len_values[k]           = len(self.variable_values[k])

        # Initialize dimensionless initial point Numpy array
        my_x0   = np.zeros([sum(len_values)])
        i       = 0
        for k in range(self.N_variables):
            prov                     = np.array(self.variable_values[k], dtype='float')
            my_x0[i:(i+len(prov))]   = prov
            i += len(prov)

        # Set bounds for DVs
        # TODO: bounds like in blade_match.py
        h = 0.4
        bds = [None]*len(my_x0)
        for i in range(len(my_x0)):
            bds[i] = (1-h, 1+h)

        # Contrainted or unconstrained optimization options
        if (self.IN['OPT_CONSTRAIN']=='NONE'):
            R = fmin_slsqp(func= self.get_objective_function, x0=my_x0, fprime=self.get_sens_objective, bounds=None, acc=1e-14, iter=20, disp=True)
        else:
            R = fmin_slsqp(func= self.get_objective_function, x0=my_x0, f_ieqcons=self.get_constrain,
                           fprime_ieqcons=self.get_sens_constrain, fprime=self.get_sens_objective, bounds=None,  acc=1e-14, iter=20, iprint=2,
                           disp=True)
        self.OPT_HIST_FILE.close()

        # Send Slack notification
        if self.SLACK_NOTIF == True:
            improv  = (R - self.RefOBJ) / self.RefOBJ * 100
            text    = "Optimisation finished after %02i steps!\nObjective function improvement [percent] = %f" % (self.DesignStep, abs(improv))
            SendSlackNotification(text)

    def run_disc_adj_validation(self):
        """
        Function to run the discrete adjoint and calculate the objective function sensitivity with respect
        to the design variables.
        """
        # Set relaxation factor to 1
        self.IN['OPT_RELAX_FAC'] = [1]

        # Create directories
        DISC_ADJ_FOLDER = 'DISCRETE_ADJOINT_'+self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']

        if self.SKIP_ADJOINT == 'NO':
            try:
                os.system("rm -rf %s"%DISC_ADJ_FOLDER)
            except:
                pass
            os.mkdir(DISC_ADJ_FOLDER)

            # Copy required files
            shutil.copyfile(DIR+self.SU2_CONFIG_IN['MESH_FILENAME'],DISC_ADJ_FOLDER+'/'+self.SU2_CONFIG_IN['MESH_FILENAME'])
            SU2_Config_change(DIR + self.SU2_CONFIG,
                            DISC_ADJ_FOLDER +'/' +
                            self.CFD_CONFIG,[],[])
            SU2_Config_change(DIR + self.SU2_CONFIG,
                            DISC_ADJ_FOLDER +'/' +
                            self.ADJ_CONFIG,['MATH_PROBLEM','CONV_FILENAME'],['DISCRETE_ADJOINT','history_ADJ'])
            shutil.copyfile(DIR + self.IN['BLADE_PAR_FILE'], DISC_ADJ_FOLDER + '/' + self.CAD_CONFIG)

            # Run tools
            os.chdir(DISC_ADJ_FOLDER)
            self.run_su2_cfd()
            self.run_adjoint()

        os.chdir(DIR + DISC_ADJ_FOLDER)
        self.run_cad_sens(self.CAD_CONFIG_IN)

        # Compute the gradient
        sens_adjoint = self.dot_product(RF=1, loc='output/sensitivities/')

        # Write output
        Fout = open('of_grad_adj.dat','w')
        Fout.write("VARIABLES=\"VARIABLE\"\t,\"GRADIENT\"\t,\"STEP\"\n")
        for i in range(len(sens_adjoint)):
            Fout.write("%i,\t%.10f,\t%.10f\n"%(i,sens_adjoint[i],0.001))
        Fout.close()

    def run_adjoint(self):
        """
        Runs SU2_CFD_AD and SU2_DOT_AD
        """
        command_AD = self.MPI_COMMAND%(self.N_CORES,'SU2_CFD_AD',self.ADJ_CONFIG)
        self.run_command(command_AD)
        command_DOT_AD = self.MPI_COMMAND%( self.N_CORES,'SU2_DOT_AD', self.ADJ_CONFIG)
        self.run_command(command_DOT_AD)


    def run_finite_diff_validation(self):
        """
        Function to obtain the objective function gradients by finite-differences
        """
        # Read user input
        try:
            FD_Step = float(IN['FIND_DIFF_STEP'][0])
        except:
            raise Exception("Finite-Differences step (FIND_DIFF_STEP) not defined in CFG file!\nExiting...")

        try:
            FD_Order = int(IN['FIND_DIFF_ORDER'][0])
        except:
            raise Exception("Finite-Differences order (FIND_DIFF_ORDER) not defined in CFG file!\nExiting...")

        if self.SKIP_REF_CFD == 'NO':
            # Check if FIND_DIFF folder exsts
            if os.path.isdir(self.FINDIFF_FOLDER):
                os.system('rm -rf FIND_DIFF')
                os.mkdir(self.FINDIFF_FOLDER)
            else:
                os.mkdir(self.FINDIFF_FOLDER)

            # Direct (reference) blade CFD is only needed for FD order equal to 1
            if FD_Order == 1:
                # Create directory
                os.mkdir(self.DIRECT_FOLDER_FD)

                # Copy required files
                shutil.copyfile(DIR + self.SU2_CONFIG,
                                self.DIRECT_FOLDER_FD + self.CFD_CONFIG)
                shutil.copyfile(DIR + self.SU2_CONFIG_IN['MESH_FILENAME'],
                                self.DIRECT_FOLDER_FD + self.SU2_CONFIG_IN['MESH_FILENAME'])

                # Run SU2_CFD for baseline blade
                os.chdir(self.DIRECT_FOLDER_FD)
                print("\n ... running SU2_CFD for reference blade ...")
                self.run_su2_cfd()

            # Create directory
            os.mkdir(self.DEFORM_FOLDER_FD)

        else:
            if os.path.isdir(self.DEFORM_FOLDER_FD):
                os.system('rm -rf FIND_DIFF/DEFORM/')
                os.mkdir(self.DEFORM_FOLDER_FD)
            else:
                os.mkdir(self.FINDIFF_FOLDER)
                os.mkdir(self.DEFORM_FOLDER_FD)

        if FD_Order == 1:
            # Read history
            if self.SKIP_REF_CFD == 'NO':
                os.chdir(self.DIRECT_FOLDER_FD)
            else:
                DISC_ADJ_FOLDER = 'DISCRETE_ADJOINT_' + self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']
                os.chdir(DISC_ADJ_FOLDER)
            ref_values = ReadHistory(self.HIST_FILENAME,1)

        # Write output file
        FD_File = open(self.FINDIFF_FOLDER + "of_grad_fd.dat", 'w')
        FD_File.write("VARIABLES=\t \"VARIABLE\",\t \"OF_GRADIENT\",\t \"FD_STEP\",\t \"SensTotalPressureLoss_1\",\t \"SensKineticEnergyLoss_1\",\t \"SensEntropyGen_1\"," \
                + "\t \"SensEulerianWork_1\",\t \"SensPressureRatio_1\",\t \"SensFlowAngleIn_1\",\t \"SensFlowAngleOut_1\",\t \"SensAbsFlowAngleIn_1\"," \
                + "\t \"SensAbsFlowAngleOut_1\",\t \"SensMassFlowIn_1\",\t \"SensMassFlowOut_1\",\t \"SensMachIn_1\",\t \"SensMachOut_1\",\t \"SensTotalEfficiency_1\"," \
                + "\t \"SensTotalStaticEfficiency_1\"\n")

        # Start loop for gradient calculation
        k = 0
        os.chdir(self.DEFORM_FOLDER_FD)
        for name in self.variable_names:

            # Length of the DV (as some variable_names might have more than one value)
            length_DV = len(self.CAD_CONFIG_IN[name])

            for j in range(length_DV):
                # Perturb design variables
                if FD_Order == 1:
                    local_hash_0    = deepcopy(self.CAD_CONFIG_IN)
                    h               = FD_Step/100*abs(local_hash_0[name][j])
                    if local_hash_0[name][j] == 0:
                        h = FD_Step/100
                    local_hash_0[name][j] += h # x+h
                    local_hash             = [local_hash_0]
                elif FD_Order == 2:
                    local_hash_0    = deepcopy(self.CAD_CONFIG_IN)
                    local_hash_1    = deepcopy(self.CAD_CONFIG_IN)
                    h               = FD_Step/100*abs(local_hash_0[name][j])
                    if local_hash_0[name][j] == 0:
                        h = FD_Step/100
                    local_hash_0[name][j] += h # x+h
                    local_hash_1[name][j] -= h # x-h
                    local_hash             = [local_hash_0, local_hash_1]
                else:
                    raise Exception("Finite-difference order not supported (FIND_DIFF_ORDER in CFG file can only be 1 or 2)!\nExiting...")

                # Initialise variable
                def_values = [None]*len(local_hash)

                for jj in range(len(local_hash)):
                    # Run write_move_surface to obtain grid dispacement
                    print("\n... starting surface calculation for variable %s run %i ..." % (name, jj))
                    self.write_move_surface(local_hash[jj])

                    # Move config file for SU2_DEF
                    shutil.copyfile(DIR + self.SU2_CONFIG,
                                    self.DEFORM_FOLDER_FD + self.DEF_CONFIG)

                    # Move initial mesh from main folder
                    shutil.copyfile(DIR + self.SU2_CONFIG_IN['MESH_FILENAME'],
                                    self.DEFORM_FOLDER_FD + self.SU2_CONFIG_IN['MESH_FILENAME'])

                    # Run SU2_DEF
                    print("\n ... running SU_DEF for variable %s run %i..." % (name, jj))
                    self.run_su2_def()

                    # Rename output .su2 mesh accordingly
                    os.rename(self.SU2_CONFIG_IN['MESH_OUT_FILENAME'], self.SU2_CONFIG_IN['MESH_FILENAME'])

                    # Copy required files for running SU2 CFD
                    shutil.copyfile(DIR + self.SU2_CONFIG,
                                    self.DEFORM_FOLDER_FD + self.CFD_CONFIG)

                    # Run SU2 CFD
                    print("\n ... running SU2_CFD for variable %s run %i..." % (name, jj))
                    self.run_su2_cfd()

                    # Gather results from CFD
                    def_values[jj] = ReadHistory(self.HIST_FILENAME,1)

                if FD_Order == 1:
                    # Forward scheme for first order FD
                    # f' = [f(x+h) - f(x)]/h
                    grads = (def_values[0] - ref_values)/h
                else:
                    # Central scheme for second order FD
                    # f' = [f(x+h) - f(x-h)]/(2h)
                    grads = (def_values[0] - def_values[1])/(2*h)

                # Find position and value of objective function
                obj_fun = grads[self.INDEX[self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']] - 1]

                # Store values in .dat file
                FD_File.write("%02i, \t%f, \t%f" % (k, obj_fun, FD_Step))
                for i in range(len(grads)):
                    FD_File.write(", \t%.10f" % grads[i])
                FD_File.write("\n")
                FD_File.flush()

                # Rename history.dat for checks
                try:
                    os.rename(self.HIST_FILENAME+".dat", self.HIST_FILENAME+ "_%s_%i.dat" % (name, k))
                except:
                    for i in range(self.NZONE):
                        hist_name = self.HIST_FILENAME + "_%i.dat"%(i)
                        os.rename(hist_name, hist_name[:-4] + "_%i.dat" % k)

                # Save volumetric grid if requested for sanity checks (only implemented for single zone)
                if self.SAVE_VOL_GRIDS and self.NZONE == 1:
                    os.rename('volumetric_grid.dat', 'volumetric_grid_%s_%i.dat' % (name, k))

                k += 1

        FD_File.close()

    def run_su2_cfd(self):
        """
        Runs SU2_CFD and SU2_SOL
        """
        command_CFD = self.MPI_COMMAND%(self.N_CORES, 'SU2_CFD', self.CFD_CONFIG)
        self.run_command(command_CFD)
        command_SOL = self.MPI_COMMAND%( self.N_CORES, 'SU2_SOL', self.CFD_CONFIG)
        self.run_command(command_SOL)

    def run_cad_sens(self, local_hash):
        """
        Runs CAD tool in SENSITIVITY mode
        :param local_hash: local hash table.
        """
        # Update DVs
        self.BladeCAD.update_DVs_control_points(deepcopy(local_hash))

        # Obtain surface coordinates and sensitivities
        self.BladeCAD.make_blade()

        # Output values
        my_plots = BladeOutput(self.BladeCAD)
        my_plots.print_blade_coordinates(path=os.getcwd())
        my_plots.print_sensitivity(path=os.getcwd())

    def run_su2_def(self):
        """
        Runs SU2_DEF
        """
        command_DEF = self.MPI_COMMAND%( self.N_CORES,'SU2_DEF', self.DEF_CONFIG)
        self.run_command(command_DEF)
    
    def write_move_surface(self,local_hash):
        """
        Writes MoveSurface.txt by comparing the surface coordinates of the deformed geometry with the ref. blade
        :param local_hash: geometry definition array
        """
        # Update values of the DVs
        self.BladeCAD.update_DVs_control_points(deepcopy(local_hash))

        # Store surface coordinates
        surface_coord_def = self.BladeCAD.get_surface_coordinates(self.u, self.v)
        surface_coord_def = np.transpose(surface_coord_def)

        # Write MoveSurface.txt
        MoveSurfOut = open('MoveSurface.txt','w')
        for i in range(len(self.match_blade)):
            diff = self.match_blade[i, 1:4] + ((surface_coord_def[i, 0:3]) - (self.surface_coord_ref[i, 0:3]))
            # pdb.set_trace()
            if self.NDIM == 2:
                MoveSurfOut.write('%i\t%.20f\t%.20f\n' % (self.match_blade[i, 0], float(diff[0]), float(diff[1])))
            else:
                MoveSurfOut.write('%i\t%.20f\t%.20f\t%.20f\n'%(self.match_blade[i,0],float(diff[0]),float(diff[1]),float(diff[2])))
        if self.DEBUG:
            if self.NDIM == 3:
                # 3D plot
                fig = plt.figure()
                ax = Axes3D(fig)
                ax.scatter(surface_coord_def[:, 0], surface_coord_def[:, 1], surface_coord_def[:, 2])
                ax.scatter(self.surface_coord_ref[:, 0], self.surface_coord_ref[:, 1], self.surface_coord_ref[:, 2])
                ax.scatter(self.match_blade[:, 1], self.match_blade[:, 2], self.match_blade[:, 3])
                ax.legend(['Deformed', 'Reference', 'Matched'])
                plt.title('3D plot')

                # 2D plot
                fig_2 = plt.figure()
                plt.plot(surface_coord_def[:, 1], surface_coord_def[:, 2], '.k', label='Deformed')
                plt.plot(self.surface_coord_ref[:, 1], self.surface_coord_ref[:, 2], '.y', label='Reference')
                plt.plot(self.match_blade[:, 2], self.match_blade[:, 3], '.r', label='Matched')
                plt.legend()
                plt.show()
            else:
                plt.plot(surface_coord_def[:, 0], surface_coord_def[:, 1], '.k', label='Deformed')
                plt.plot(self.surface_coord_ref[:, 0], self.surface_coord_ref[:, 1], '.y', label='Reference')
                plt.plot(self.match_blade[:, 1], self.match_blade[:, 2], '.r', label='Matched')
                plt.legend()
                plt.axis('equal')
                plt.show()

        MoveSurfOut.close()

    def store_ref_blade(self):
        """
        Function that stores in global variables the baseline blade properties
        """
        print("... starting calculations for reference blade ...")
        try:
            self.match_blade = np.loadtxt(DIR + self.MATCH_BLADE, delimiter=',', skiprows=1)
        except:
            raise Exception('Match blade filename in CFG file does not exist in the directory... or error reading file...\nExiting...')

        # Create baseline blade class
        self.BladeBase = Blade3D(deepcopy(self.CAD_CONFIG_IN))

        # Update u and v
        self.u = np.where(self.match_blade[:, -2] > 1, self.match_blade[:, -2] - 1, self.match_blade[:, -2])
        self.v = self.match_blade[:, -1]
        self.BladeBase.update_uv_values(self.u, self.v)

        # Calculate surface coordinates and sensitivities
        # self.BladeBase.make_surface_interpolant(interp_method='bilinear')
        self.surface_coord_ref = self.BladeBase.get_surface_coordinates(self.u, self.v)
        self.surface_coord_ref = np.transpose(self.surface_coord_ref)

        # Create blade class for RunCAD blades
        self.BladeCAD   = Blade3D(deepcopy(self.CAD_CONFIG_IN))
        self.BladeCAD.update_uv_values(self.u, self.v)

    def get_sens_objective(self,parameter):
        """
        Calculates the sensitivity of the objective function
        """
        if self.RESTART == 'YES' and os.path.exists(self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']) + 'sens_obj.pickle'):
            # Change directory to discrete adjoint folder
            os.chdir(self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']))

            # Load sensitivity from the pickle file (sensitivity is already multiplied by RF)
            self.OBJ_SENS_PCKL = open("sens_obj.pickle", "rb")
            total_sensitivity = pickle.load(self.OBJ_SENS_PCKL)
            if self.DesignStep == 1:
                self.RefSens = total_sensitivity
            print("Sensitivity (old) for design step no. %02i: " % self.DesignStep)
            print(total_sensitivity)

            # Multiply sensitivity by relaxation factor (this is required in case there is a change of RF between steps)
            if (self.DesignStep is not 1) and (self.DesignStep is not self.DesignStep_Restart):
                total_sensitivity = total_sensitivity / self.of_relax_factor[self.DesignStep - 1] * self.of_relax_factor[self.DesignStep]
            elif self.DesignStep is self.DesignStep_Restart:
                # If the current design step coincides with the last available design step
                # multiply the sensitivity by the relaxation factor set in the configuration file
                total_sensitivity = total_sensitivity / self.of_relax_factor[self.DesignStep - 1] * float(self.IN['OPT_RELAX_FAC'][0])

            print("Sensitivity (new) for design step no. %02i: " % self.DesignStep)
            print(total_sensitivity)

        else:
            # Create directories
            os.mkdir(self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']))
            os.chdir(self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']))

            # Run SU2 discrete adjoint
            self.Module = "ADJ_OBJ"
            self.move_files()
            if os.path.exists(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/' + self.SENS_FILENAME) and self.DesignStep == 1:
                # Copy history_ADJ.dat and flow.dat if a previous adjoint run exists
                shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/' + self.SENS_FILENAME,
                                self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']) + '/' + self.SENS_FILENAME)
                if self.NZONE == 1:
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/history_ADJ.dat',
                                    self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']) + '/history_ADJ.dat')
                elif self.NZONE == 2:
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/history_ADJ_0.dat',
                                    self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN[
                                        'OBJECTIVE_FUNCTION']) + '/history_ADJ_0.dat')
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/history_ADJ_1.dat',
                                    self.SENS_OBJ_FOLDER % (self.DesignStep, self.SU2_CONFIG_IN[
                                        'OBJECTIVE_FUNCTION']) + '/history_ADJ_1.dat')

            else:
                self.run_adjoint()

            # Perform dot product between CAD and adjoint sensitivities
            total_sensitivity = self.dot_product(RF=float(self.IN['OPT_RELAX_FAC'][0]))

            if self.DesignStep == 1:
                self.RefSens = total_sensitivity

            # Write sensitivities (with relax. factor) in the pickle file
            self.OBJ_SENS_PCKL = open("sens_obj.pickle", "wb")
            pickle.dump(total_sensitivity, self.OBJ_SENS_PCKL)
            self.OBJ_SENS_PCKL.close()

        print(total_sensitivity)

        return total_sensitivity

    def get_sens_constrain(self,parameter):
        """
        Calculates the sensitivity of the constraint (OPT_CONSTRAIN) defined in the input file
        """
        if self.RESTART == 'YES' and os.path.exists(self.SENS_OBJ_FOLDER%(self.DesignStep,self.IN['OPT_CONSTRAIN'][0]) + 'sens_cons.pickle'):
            # Change directory
            os.chdir(self.SENS_OBJ_FOLDER % (self.DesignStep, self.IN['OPT_CONSTRAIN'][0]))

            # Load sensitivity from the pickle file (sensitivity is already multiplied by RF)
            self.CONS_SENS_PCKL = open("sens_cons.pickle","rb")
            total_sensitivity = pickle.load(self.CONS_SENS_PCKL)

            # Multiply sensitivity by relaxation factor (required if different RF have been used)
            if (self.DesignStep is not 1) and (self.DesignStep is not self.DesignStep_Restart):
                total_sensitivity = total_sensitivity / self.cons_relax_factor[self.DesignStep - 1] * self.cons_relax_factor[self.DesignStep]
            elif self.DesignStep is self.DesignStep_Restart:
                # If the current design step coincides with the last available design step
                # multiply the sensitivity by the relaxation factor set in the configuration file
                total_sensitivity = total_sensitivity / self.cons_relax_factor[self.DesignStep - 1] * float(self.IN['CONS_RELAX_FAC'][0])

        else:
            # Create directories
            os.mkdir(self.SENS_OBJ_FOLDER%(self.DesignStep,self.IN['OPT_CONSTRAIN'][0]))
            os.chdir(self.SENS_OBJ_FOLDER%(self.DesignStep,self.IN['OPT_CONSTRAIN'][0]))

            # Run SU2 discrete adjoint
            if self.DesignStep == 1 and os.path.exists(DIR + 'DISCRETE_ADJOINT_FLOW_ANGLE_OUT/'):
                shutil.copyfile(DIR + 'DISCRETE_ADJOINT_FLOW_ANGLE_OUT/surface_sens.dat', self.SENS_OBJ_FOLDER % (
                self.DesignStep, self.IN['OPT_CONSTRAIN'][0]) + 'surface_sens.dat')
                os.chdir(self.SENS_OBJ_FOLDER % (self.DesignStep, self.IN['OPT_CONSTRAIN'][0]))
            else:
                self.Module = "ADJ_CON"
                self.move_files()
                self.run_adjoint()
            total_sensitivity = self.dot_product(RF=float(self.IN['CONS_RELAX_FAC'][0]))

            # Write sensitivity in pickle file
            self.CONS_SENS_PCKL = open("sens_cons.pickle", "wb")
            pickle.dump(total_sensitivity, self.CONS_SENS_PCKL)
            self.CONS_SENS_PCKL.close()

        print(total_sensitivity)

        return np.array([total_sensitivity])


    def dot_product(self,RF=1,loc='../CAD/output/sensitivities/'):
        """
        Performs DOT product between CAD and CFD and returns sensitivity with respect to the Design Variables
        :param loc: Location of the CAD sensitivities files
        :param RF: Relaxation factor
        :return: Array with sensitivities with respect to the design variables
        """
        # Read CFD sensitivities
        len_sens        = file_length(self.SENS_FILENAME)
        stop_sens       = file_endread(self.SENS_FILENAME)
        sensitivity_CFD = np.genfromtxt(self.SENS_FILENAME, dtype='float', skip_header=3, skip_footer=len_sens-stop_sens)

        # Define positions of interest in the CFD sensitivities
        x_sens_pos = int(self.NDIM) + 0
        y_sens_pos = int(self.NDIM) + 1
        z_sens_pos = int(self.NDIM) + 2

        total_sensitivity  = [None] * self.N_DVs
        sens_index         = 0

        # Loop for all design variables
        for name in self.variable_names:

            # Read CAD sensitivities
            length_DV = len(self.CAD_CONFIG_IN[name])
            for k in range(length_DV):
                # Initialize and read sensitivities
                sensitivity_var = 0
                sensitivity_CAD = np.genfromtxt(loc+"/grad_"+name+"_%i.csv" % k,skip_header=1, delimiter=',')
                if self.DEBUG and k == 0:
                    if self.NDIM == 2:
                        plt.plot(sensitivity_CAD[:, 1], sensitivity_CAD[:, 2], '.r')
                        plt.plot(sensitivity_CFD[:, 0], sensitivity_CFD[:, 1], '.k')
                        plt.axis('equal')
                        plt.legend(['CAD matched', 'CFD'])
                        plt.show()
                    else:
                        # 3D plot
                        fig = plt.figure()
                        ax = Axes3D(fig)
                        ax.scatter(sensitivity_CAD[:,3], sensitivity_CAD[:,2], sensitivity_CAD[:,1])
                        ax.scatter(sensitivity_CFD[:,2], sensitivity_CFD[:,1], sensitivity_CFD[:,0])
                        ax.set_xlabel('Z-axis')
                        ax.set_ylabel('Y-axis')
                        ax.set_zlabel('X-axis')
                        ax.legend(['CAD matched', 'CFD'])
                        plt.title('3D plot')
                        # 2D plane
                        fig2 = plt.figure()
                        plt.plot(sensitivity_CAD[:,3], sensitivity_CAD[:,2], '.')
                        plt.plot(sensitivity_CFD[:, 2], sensitivity_CFD[:, 1], '.')
                        plt.legend(['CAD matched', 'CFD'])
                        plt.grid()
                        plt.title('Z-Y plane')
                        plt.show()
                for i in range(len(sensitivity_CFD)):
                    # Calculate difference in X, Y and Z (3D cases) coord.
                    diff_x = sensitivity_CAD[:,1] - sensitivity_CFD[i,0]
                    diff_y = sensitivity_CAD[:,2] - sensitivity_CFD[i,1]
                    if self.NDIM == 2:
                        diff_z = 0.0
                    else:
                        diff_z = sensitivity_CAD[:,3] - sensitivity_CFD[i,2]

                    # Calculate distance from CFD to CAD point
                    dist = (diff_x ** 2 + diff_y ** 2 + diff_z ** 2) ** 0.5

                    # Find position of minimum distance
                    index_min = np.argmin(dist)
                    # print(dist[index_min])

                    # Compute the sensitivities X, Y and Z (3D cases)
                    sensitivity_x = sensitivity_CAD[index_min, x_sens_pos + 1] * sensitivity_CFD[i, x_sens_pos]            # x-total_sensitivity
                    sensitivity_y = sensitivity_CAD[index_min, y_sens_pos + 1] * sensitivity_CFD[i, y_sens_pos]            # y-total_sensitivity
                    if self.NDIM == 2:
                        sensitivity_z = 0.0
                    else:
                        sensitivity_z = sensitivity_CAD[index_min, z_sens_pos + 1] * sensitivity_CFD[i, z_sens_pos]        # z-total_sensitivity

                    sensitivity_tot = (sensitivity_x ** 2 + sensitivity_y ** 2 + sensitivity_z ** 2) ** 0.5                                     # total
                    sensitivity_var += sensitivity_x + sensitivity_y + sensitivity_z                                                     # dot product

                # Multiply sensitivities by relaxation factor (required for optimisation)
                total_sensitivity[sens_index] = sensitivity_var * RF
                sens_index += 1

        print("... SENSITIVITY SUMMARY (after dot product) ...")
        print(total_sensitivity)

        # Store reference sensitivities array
        if self.DesignStep == 1:
            self.RefSens = np.array(total_sensitivity)

        return np.array(total_sensitivity)

    def make_local_hash(self,parameter):
        """
        Writes local hash necessary to run the CAD tool for perturbed DVs
        :param parameter: design variables
        :return: local hash table
        """
        local_index = 0
        local_hash = deepcopy(self.CAD_CONFIG_IN)
        for i in self.variable_names:
            for k in range(len(self.CAD_CONFIG_IN[i])):
                local_hash[i][k] = parameter[local_index]
                local_index += 1

        # Output local hash to a config file
        self.write_parablade_cfg(local_hash)

        return local_hash

    def get_objective_function(self,parameter):
        """
        Function to obtain the objective function value. Runs CAD geometry tool, SU2_DEF and SU2_CFD
        :return: Dimensionless objective function value
        """
        # Used by the Optimizer Only
        # Create directories
        self.DesignStep += 1

        # Change to 'False' the REMESH boolean
        if self.REMESH:
            if self.DesignStep > self.REMESH_STEP:
                self.REMESH = False
        #pdb.set_trace()
        if self.RESTART == 'YES' and os.path.exists(self.DIRECT_FOLDER % self.DesignStep + 'obj_fun.pickle'):
            # Change directory
            os.chdir(self.DIRECT_FOLDER % self.DesignStep)

            # Load pickle file containing objective function value
            self.OBJ_PCKL = open("obj_fun.pickle", "rb")
            obj_fun = pickle.load(self.OBJ_PCKL)
            if self.DesignStep == 1:
                self.RefOBJ = obj_fun

            # Display optset pickle to check if it's the same (useful for debugging)
            self.OPTSET_PCKL = open(self.DESIGN_FOLDER % self.DesignStep + "optset.pickle", 'rb')
            var_nams = pickle.load(self.OPTSET_PCKL)
            var_vals = pickle.load(self.OPTSET_PCKL)
            par = pickle.load(self.OPTSET_PCKL)
            print("... correspondance between DVs (restarted vs. decided by optimizer) ...")
            print(par == parameter)

        else:
            # Used by the Optimizer Only
            # Create directories
            self.make_design_folder()

            # Create and write pickle file for possible restoring containing optimization settings
            self.OPTSET_PCKL = open(self.DESIGN_FOLDER % self.DesignStep + "optset.pickle", "wb")
            pickle.dump(self.variable_names, self.OPTSET_PCKL)                      # Variable names
            pickle.dump(self.variable_values, self.OPTSET_PCKL)                     # Variable values
            pickle.dump(parameter, self.OPTSET_PCKL)                                # Design vector
            pickle.dump(float(self.IN['OPT_RELAX_FAC'][0]), self.OPTSET_PCKL)       # Optimization relax. factor
            if self.IN['OPT_CONSTRAIN'] is not 'NONE':
                pickle.dump(float(self.IN['CONS_RELAX_FAC'][0]), self.OPTSET_PCKL)  # Constrain relax. factor
            self.OPTSET_PCKL.close()

            # ONLY FOR REMESHING
            # Will read new matched blade and update U and V accordingly for the base blade and CAD blade
            if self.REMESH and self.DesignStep == self.REMESH_STEP:
                print("... updating u- and v- for the new remeshed blade ...")
                # Read new matched file
                self.match_blade = np.loadtxt(DIR + '/EXTRA/REMESH_1/' + self.MATCH_BLADE, delimiter=',', skiprows=1)
                # Sort u and v
                self.u = np.where(self.match_blade[:, -2] > 1, self.match_blade[:, -2] - 1,
                                  self.match_blade[:, -2])
                self.v = self.match_blade[:, -1]
                # Update u and v for the reference blade
                self.BladeBase.update_uv(self.u, self.v)
                # Calculate surface coordinates for the new u and v values (ref blade)
                self.surface_coord_ref = self.BladeBase.get_surface_coordinates(self.u, self.v)
                # self.surface_coord_ref = np.transpose(self.surface_coord_ref)
                # Update u and v for the CAD blade
                self.BladeCAD.update_uv(self.u, self.v)
                print("... finished updating u- and v- for remeshed blade ...")

            # Define geometry
            local_hash = self.make_local_hash(parameter)

            # Obtain mesh displacement only for DesignStep != 1
            os.chdir(self.CAD_FOLDER % self.DesignStep)
            if self.DesignStep == 1:
                pass
            else:
                self.write_move_surface(local_hash)

            # Run CAD to obtain geometrical sensitivities
            self.run_cad_sens(local_hash)

            # Run SU2_DEF only for DesignStep != 1
            # For DesignStep = 1 SU2_DEF will be run to generate a MoveSurface.txt for the similar mesh search
            # (only if MESH_SEARCH is set to TRUE)
            if self.DesignStep == 1:
                if self.MESH_SEARCH:
                    os.chdir(self.DEFORM_FOLDER % self.DesignStep)
                    # Move config file
                    shutil.copyfile(DIR + self.SU2_CONFIG,
                                    self.DEFORM_FOLDER % self.DesignStep + self.DEF_CONFIG)
                    # Move mesh file
                    shutil.copyfile(DIR + self.SU2_CONFIG_IN['MESH_FILENAME'],
                                    self.DEFORM_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME'])
                    # Run SU2_DEF to generate MoveSurface.txt
                    # (run on single core as multicore can generate corrupt MoveSurface.txt files)
                    command_DEF = "SU2_DEF " + self.DEF_CONFIG
                    self.run_command(command_DEF)
                else:
                    pass
            else:
                os.chdir(self.DEFORM_FOLDER % self.DesignStep)
                self.Module = "DEF"
                self.move_files()
                if not self.REMESH:
                    self.run_su2_def()
                elif self.REMESH and self.DesignStep == self.REMESH_STEP:
                    os.rename(self.SU2_CONFIG_IN['MESH_FILENAME'], self.SU2_CONFIG_IN['MESH_OUT_FILENAME'])

            # Run SU2_CFD
            os.chdir(self.DIRECT_FOLDER % self.DesignStep)
            self.Module = "CFD"
            self.move_files()
            if os.path.exists(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/' + self.HIST_FILENAME + '.dat') and self.DesignStep == 1:
                # Copy history.dat, flow.dat and restart_flow.dat if a previous adjoint run exists
                shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/' + self.HIST_FILENAME + '.dat',
                                self.DIRECT_FOLDER % self.DesignStep + '/' + self.HIST_FILENAME + '.dat')
                if self.NZONE == 1:
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/flow.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/flow.dat')
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/restart_flow.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/restart_flow.dat')
                elif self.NZONE == 2:
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/flow_0.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/flow_0.dat')
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/flow_1.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/flow_1.dat')
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/restart_flow_0.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/restart_flow_0.dat')
                    shutil.copyfile(DIR + 'DISCRETE_ADJOINT_' + self.OBJ_FUNCTION + '/restart_flow_1.dat',
                                    self.DIRECT_FOLDER % self.DesignStep + '/restart_flow_1.dat')
            else:
                self.run_su2_cfd()

            # Read converged solution and objective function value
            conv    = ReadHistory(self.HIST_FILENAME,self.NZONE)
            obj_fun = conv[self.INDEX[self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']] - 1]

            # Store initial objective function value
            if self.DesignStep == 1:
                self.RefOBJ = obj_fun

            # Obtain constrain value
            # For unconstrained problems this will be equal to 0.0
            print("Called in get_objective_function ...")
            constr = self.get_constrain(parameter)

            # Write optimization history file
            self.OPT_HIST_FILE.write("%i\t,%.10f\t,%.10f\n" % (self.DesignStep, obj_fun, constr))
            self.OPT_HIST_FILE.flush()

            # Write objective function to pickle file
            self.OBJ_PCKL = open("obj_fun.pickle", "wb")
            pickle.dump(obj_fun, self.OBJ_PCKL)
            self.OBJ_PCKL.close()

        # Output to Slack
        if self.SLACK_NOTIF == True:
            improv  = (obj_fun - self.RefOBJ)/self.RefOBJ*100
            text    = "Design step: %02i       Objective function improvement [percent] = %f" % (self.DesignStep, abs(improv))
            SendSlackNotification(text)

        print('... Objective function value (non-dimensional) ...')
        print(obj_fun/self.RefOBJ)

        return obj_fun/self.RefOBJ

    def get_constrain(self, parameter):
        """
        Function to get the constraint value
        :return: Constraint value (except for unconstrained optimizations, where constr = 0.0)
        """

        if self.RESTART == 'YES' and os.path.exists(self.DIRECT_FOLDER % (self.DesignStep) + 'cons.pickle'):
            try:
                os.chdir(self.DIRECT_FOLDER % (self.DesignStep))
                self.CONS_PCKL = open("cons.pickle", "rb")
                constr = pickle.load(self.CONS_PCKL)
            except:
                return [0]

        else:
            # Used by the Optimizer Only
            try:
                os.chdir(self.DIRECT_FOLDER % (self.DesignStep))
            except:
                return([0])
            history = np.loadtxt(self.HIST_FILENAME + '.dat', skiprows=3, delimiter=',')

            if self.IN['OPT_CONSTRAIN'][1] is '>':
                constr_abs = abs(history[-1][self.INDEX[self.IN['OPT_CONSTRAIN'][0]]])
                constr = constr_abs - self.IN['OPT_CONSTRAIN'][2]
            elif self.IN['OPT_CONSTRAIN'][1] is '<':
                constr_abs = abs(history[-1][self.INDEX[self.IN['OPT_CONSTRAIN'][0]]])
                constr = self.IN['OPT_CONSTRAIN'][2] - constr_abs
            elif self.IN['OPT_CONSTRAIN'][1] is ':':
                h = 0.02  # Maximum allowed deviation from the equality constrain
                constr_abs = abs(history[-1][self.INDEX[self.IN['OPT_CONSTRAIN'][0]]])
                constr_ub = 1 - constr_abs / (self.IN['OPT_CONSTRAIN'][2] * (1 + h))  # Upper bound
                constr_lb = constr_abs / (self.IN['OPT_CONSTRAIN'][2] * (1 - h)) - 1  # Lower bound
                constr = np.array([constr_ub, constr_lb])
            elif self.IN['OPT_CONSTRAIN'][1] is 'O':  # If OPT_CONSTRAIN = NONE
                constr_abs = abs(history[-1][self.INDEX['FLOW_ANGLE_OUT']])
                constr = constr_abs
            else:
                raise Exception("OPT_CONSTRAIN in Optimization.cfg wrongly defined!\nExiting...")

            # Write constrain in pickle file
            self.CONS_PCKL = open("cons.pickle", "wb")
            pickle.dump(constr, self.CONS_PCKL)
            self.CONS_PCKL.close()

        if self.DesignStep == 1:
            self.RefCons = constr

        # Output to Slack
        if self.SLACK_NOTIF == True:
            text = "Design step: %02i       Constraint value = %f" % (self.DesignStep, constr)
            SendSlackNotification(text)

        print('... Constraint value ...')
        print(constr)

        return np.array([constr])

    def make_design_folder(self):
        """
        Function to create the folder directory for the Shape Optimization case
        """
        os.mkdir(self.DESIGN_FOLDER%(self.DesignStep))
        os.mkdir(self.DEFORM_FOLDER%(self.DesignStep))
        os.mkdir(self.CAD_FOLDER%(self.DesignStep))
        os.mkdir(self.DIRECT_FOLDER%(self.DesignStep))

    def write_parablade_cfg(self, local_hash):
        """
        Function that writes the CAD configuration file (BladeOpti.cfg) for keeping track of the geometry changes
        :param local_hash: dictionary containing the geometry information
        """
        # Define name and location of the file
        parablade_cfg = self.CAD_FOLDER % (self.DesignStep) + self.CAD_CONFIG

        # Write file
        fileOUT = open(parablade_cfg, 'w')
        WriteConfigFile(fileOUT, local_hash)
        fileOUT.close()

    def move_files(self):
        """
        Function to move REQUIRED filed by each module. Firstly, the function classifies the files that need to be moved
        inside an IF statement and then it moves them outside the conditional.
        """
        if self.Module == "CFD":
            # Config files
            first_out  = DIR + self.SU2_CONFIG
            first_in   = self.DIRECT_FOLDER % self.DesignStep + self.CFD_CONFIG

            # Mesh files
            if not self.REMESH:
                if self.DesignStep == 1:
                    second_out = DIR + self.SU2_CONFIG_IN['MESH_FILENAME']
                else:
                    second_out = self.DEFORM_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_OUT_FILENAME']
            elif self.REMESH and self.DesignStep == self.REMESH_STEP:
                second_out   = DIR + '/EXTRA/REMESH_1/' + self.SU2_CONFIG_IN['MESH_FILENAME']
            elif self.REMESH and self.DesignStep < self.REMESH_STEP:
                if self.DesignStep == 1:
                    second_out = DIR + self.SU2_CONFIG_IN['MESH_FILENAME']
                else:
                    second_out   = self.DEFORM_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_OUT_FILENAME']
            second_in  = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']

            # extra_files files
            extra_files = False

            # Config change
            change_cfg = False

        elif self.Module == "DEF":
            # Config files
            first_out  = DIR + self.SU2_CONFIG
            first_in   = self.DEFORM_FOLDER % self.DesignStep + self.DEF_CONFIG

            # Mesh files
            # If MESH_SEARCH is not activated it will move the mesh from the previous design step
            if not self.REMESH:
                if not self.MESH_SEARCH and self.DesignStep > 2:
                    second_out = self.DEFORM_FOLDER % (self.DesignStep - 1) + self.SU2_CONFIG_IN['MESH_OUT_FILENAME']
                elif (self.MESH_SEARCH and self.DesignStep > 2):
                    second_out = self.search_similar_mesh()
                else:
                    second_out = DIR + self.SU2_CONFIG_IN['MESH_FILENAME']
            else:
                if self.DesignStep == self.REMESH_STEP:
                    second_out = DIR + '/EXTRA/REMESH_1/' +  \
                                self.SU2_CONFIG_IN['MESH_FILENAME']
                elif self.DesignStep < self.REMESH_STEP:
                    second_out = DIR + self.SU2_CONFIG_IN['MESH_FILENAME']
            second_in  = self.DEFORM_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']

            # extra_files files - MoveSurface
            extra_files    = True
            extra_out = self.CAD_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MOTION_FILENAME']
            extra_in  = self.DEFORM_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MOTION_FILENAME']

            # Config change
            change_cfg = False

        elif self.Module == "ADJ_OBJ":
            OBJ_FUN     = self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']

            # Flow solution files
            first_out    = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
            first_in     = self.SENS_OBJ_FOLDER % (self.DesignStep, OBJ_FUN) + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']

            # Mesh files
            if not self.REMESH:
                second_out = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']
            elif self.REMESH and self.DesignStep == self.REMESH_STEP:
                second_out = DIR + '/EXTRA/REMESH_1/' + self.SU2_CONFIG_IN['MESH_FILENAME']
            elif self.REMESH and self.DesignStep < self.REMESH_STEP:
                second_out = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']
            second_in    = self.SENS_OBJ_FOLDER % (self.DesignStep, OBJ_FUN) + self.SU2_CONFIG_IN['MESH_FILENAME']

            # extra_files files
            if self.NZONE == 1:
                extra_files   = False
            elif self.NZONE == 2:
                # Multizone files - flow solution for second zone
                extra_files  = True
                extra_out    = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
                extra_in     = self.SENS_OBJ_FOLDER % (self.DesignStep, OBJ_FUN) + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
                # Changes in filename ending for multizone
                first_out    = first_out[:-4] + "_1.dat"
                first_in     = first_in[:-4] + "_1.dat"
                extra_out    = extra_out[:-4] + "_0.dat"
                extra_in     = extra_in[:-4] + "_0.dat"

            # Config change
            change_cfg   = True
            ref_cfg      = DIR + self.SU2_CONFIG
            new_cfg      = self.SENS_OBJ_FOLDER % (self.DesignStep, OBJ_FUN) + self.ADJ_CONFIG
            if self.NZONE == 1:
                option_cfg   = ['MATH_PROBLEM', 'MARKER_DEFORM_TANGENTIAL']
                new_option_cfg   = ['DISCRETE_ADJOINT', 'NONE']
            elif self.NZONE == 2:
                # TODO: hardcoded for APU test case
                option_cfg   = ['MATH_PROBLEM', 'RAMP_OUTLET_PRESSURE', 'RAMP_ROTATING_FRAME','MARKER_DEFORM_TANGENTIAL']
                new_option_cfg   = ['DISCRETE_ADJOINT', 'NO', 'NO','NONE']

        elif self.Module == "ADJ_CON":
            CONSTR      = self.IN['OPT_CONSTRAIN'][0]

            # Flow solution files
            first_out    = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
            first_in     = self.SENS_OBJ_FOLDER % (self.DesignStep, CONSTR) + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']

            # Mesh files
            if not self.REMESH:
                second_out   = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']
            elif self.REMESH and self.DesignStep == self.REMESH_STEP:
                second_out = DIR + '/EXTRA/REMESH_1/' + self.SU2_CONFIG_IN['MESH_FILENAME']
            elif self.REMESH and self.DesignStep < self.REMESH_STEP:
                second_out   = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MESH_FILENAME']
            second_in    = self.SENS_OBJ_FOLDER % (self.DesignStep, CONSTR) + self.SU2_CONFIG_IN['MESH_FILENAME']

            # extra_files files
            if self.NZONE == 1:
                extra_files   = False
            elif self.NZONE == 2:
                # Multizone files - flow solution for second zone
                extra_files       = True
                extra_out    = self.DIRECT_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
                extra_in     = self.SENS_OBJ_FOLDER % (self.DesignStep, CONSTR) + self.SU2_CONFIG_IN['SOLUTION_FLOW_FILENAME']
                # Changes in filename ending for multizone
                first_out    = first_out[:-4] + "_1.dat"
                first_in     = first_in[:-4] + "_1.dat"
                extra_out    = extra_out[:-4] + "_0.dat"
                extra_in     = extra_in[:-4] + "_0.dat"

            # Config change
            change_cfg   = True
            ref_cfg      = DIR + self.SU2_CONFIG
            new_cfg      = self.SENS_OBJ_FOLDER % (self.DesignStep, CONSTR) + self.ADJ_CONFIG
            if self.NZONE == 1:
                option_cfg   = ['MATH_PROBLEM', 'OBJECTIVE_FUNCTION']
                new_option_cfg   = ['DISCRETE_ADJOINT', CONSTR]
            elif self.NZONE == 2:
                # TODO: hardcoded for APU test case
                option_cfg   = ['MATH_PROBLEM', 'OBJECTIVE_FUNCTION', 'RAMP_OUTLET_PRESSURE', 'RAMP_ROTATING_FRAME']
                new_option_cfg   = ['DISCRETE_ADJOINT', CONSTR, 'NO', 'NO']

        # Move first file
        shutil.copyfile(first_out, first_in)

        # Move second file
        # Creates a symbolic link for the mesh files - allows for memory saving
        os.symlink(second_out, second_in)

        # Move extra files
        if extra_files:
            shutil.copyfile(extra_out, extra_in)

        # Change and move config file
        if change_cfg:
            SU2_Config_change(ref_cfg, new_cfg, option_cfg, new_option_cfg)

    def run_plot_validation(self, loc_ADJ = 'DISCRETE_ADJOINT_' + IN['OBJECTIVE_FUNCTION'], loc_FD = FINDIFF_FOLDER):
        """
        Function to plot results of the gradient validation.
        Discrete adjoint and FD modules need to be ran beforehand.
        :param: loc_ADJ = location of the adjoint simulations
        :param: loc_FD = location of the finite difference simulations
        """
        # Read adjoint results
        try:
            adjoint_gradients = np.genfromtxt(loc_ADJ + "/of_grad_adj.dat", delimiter=',', dtype='float',
                                     skip_header=1)
        except:
            raise Exception("File of_grad_adj.dat non-existing in the folder!\nExiting...")

        try:
            findiff_gradients = np.genfromtxt(loc_FD + "/of_grad_fd.dat", delimiter=',', dtype='float',
                                     skip_header=1)
            if findiff_gradients.ndim == 1:
                n_DVs = np.array([1])
                findiff_gradients = [findiff_gradients]
            else:
                n_DVs = np.array(range(len(findiff_gradients[:, 0])))
        except:
            raise Exception("File of_grad_fd.dat non existing in the folder!\nExiting...")

        # Plot selected variable names if ALL the DVs are written in the sensitivity files
        if (len(adjoint_gradients[:,1]) == len(self.section_variable_names)\
                and len(findiff_gradients[:,1]) == len(self.section_variable_names)):
            dict_ADJ    = dict(zip(self.section_variable_names, adjoint_gradients[:,1]))
            dict_FD     = dict(zip(self.section_variable_names, findiff_gradients[:,1]))

            Grads_FD_plot = np.zeros(len(self.variable_names))
            Grads_ADJ_plot = np.zeros(len(self.variable_names))

            i = 0
            for name in self.variable_names:
                Grads_FD_plot[i]    = dict_FD[name]
                Grads_ADJ_plot[i]   = dict_ADJ[name]
                i += 1

            n_DVs = np.array(range(len(Grads_FD_plot)))

            xTicksChange = True
        else:
            xTicksChange = False
            Grads_FD_plot = findiff_gradients[:,1]
            Grads_ADJ_plot = adjoint_gradients[:,1]


        # Read FD step for title plot
        FD_Step = findiff_gradients[0][2]

        # Set 'usetex' to True to have LaTeX fonts in the figures
        mpl.rc('text', usetex=False)
        mpl.rc('font', family='serif')

        # Plots
        fig_1 = plt.figure()
        if n_DVs.size != 1:
            plt.plot(n_DVs, Grads_FD_plot, 'k*', label='FD ' + r"$\Delta h =  %s$ %%" % str(FD_Step))
            plt.plot(n_DVs, Grads_ADJ_plot, 'r^', label ='ADJ')
            if xTicksChange:
                x_ticks = self.variable_names
                plt.xticks(n_DVs, x_ticks, rotation=45)
        else:
            plt.plot(n_DVs, Grads_FD_plot, 'k*', label='FD ' + r"$\Delta h =  %s$ %%" % str(FD_Step))
            plt.plot(n_DVs, Grads_ADJ_plot, 'r^', label ='ADJ')
        # plt.title(r"FD vs AD Gradient Validation")
        plt.xlabel(r"$\alpha$ [-]", fontsize=12)
        plt.ylabel(r"$\frac{dJ}{d\alpha}$ [-]", fontsize=12)
        plt.legend()
        plt.grid(True)
        plt.savefig('grad_val_ad_fd_DV_h%s.eps' % str(FD_Step),bbox_inches='tight',format='eps', dpi=2000)


        # Different style plot
        error = 0.20             # Error line to be plotted

        fig_2 = plt.figure()
        x_min_FD = np.amin(Grads_FD_plot)
        x_max_FD = np.amax(Grads_FD_plot)
        x_range_FD = np.linspace(1.5 * x_min_FD, 1.4 * x_max_FD, 501)

        # Error line
        x_err    = np.linspace(1.5 * x_min_FD, 1.4 * x_max_FD, 501)
        x_err_2  = np.linspace(1.5 * x_min_FD, 1.4 * x_max_FD, 501)
        x_err[np.where(x_err < 0.0)] = (1-error) * x_err[np.where(x_err < 0.0)]
        x_err[np.where(x_err > 0.0)] = (1+error) * x_err[np.where(x_err > 0.0)]
        x_err_2[np.where(x_err_2 < 0.0)] = (1+error) * x_err_2[np.where(x_err_2 < 0.0)]
        x_err_2[np.where(x_err_2 > 0.0)] = (1-error) * x_err_2[np.where(x_err_2 > 0.0)]

        plt.plot(Grads_ADJ_plot, Grads_FD_plot, 'o')
        plt.plot(x_range_FD, x_range_FD, 'k--',
                 label=r"$(\frac{d J}{d\alpha})_{AD} = (\frac{d J}{d\alpha})_{FD}$")
        plt.plot(-x_range_FD, -x_range_FD, 'k--')
        plt.plot(x_range_FD, x_err, 'r--',
                 label=r"$(\frac{d J}{d\alpha})_{AD} = (\frac{d J}{d\alpha})_{FD}" + chr(177) + r"%s (\frac{d J}{d\alpha})_{FD}$" % (error))
        plt.plot(x_range_FD, x_err_2, 'r--')
        plt.plot(-x_range_FD, -x_err, 'r--')
        plt.plot(-x_range_FD, -x_err_2, 'r--')
        plt.xlabel(r"$(\frac{d J}{d\alpha})_{AD}$ [-]", fontsize=12)
        plt.ylabel(r"$(\frac{d J}{d\alpha})_{FD}$ [-]", fontsize=12)
        plt.axis('equal')
        plt.legend()
        plt.savefig('grad_val_ad_fd_h%s.eps' % str(FD_Step),bbox_inches='tight',format='eps', dpi=2000)
        if self.SHOW_PLOT:
            plt.show()

    def design_variable_checks(self):
        """
        Function that performs checks for consistency for a given set of DVs
        """
        #TODO: needs to be corrected with new names

        # In 2D only section DVs can be used
        if self.NDIM == 2:
            for elem in self.variable_names:
                if not elem in self.section_variable_names:  # Checks for consistency
                    if elem in self.channel_variable_names:
                        print("Use of channel DVs ('%s')" % elem + " for 2D cases not supported")
                    else:
                        print("Design variable '%s'" % elem + " not valid!")
                    print("Using default design variables...\nPress Ctrl+C to abort...")
                    try:
                        time.sleep(5)
                        self.variable_names = self.section_variable_names
                    except KeyboardInterrupt:
                        raise Exception("\nExiting...")

        # In 3D both channel and section DVs can be used
        elif self.NDIM == 3:
            for elem in self.variable_names:
                if not elem in (self.section_variable_names + self.channel_variable_names):  # Checks for consistency
                    print("Design variable '%s'" % elem + " not valid!")
                    print("Using default design variables...\nPress Ctrl+C to abort...")
                    try:
                        time.sleep(5)
                        self.variable_names = self.section_variable_names + self.channel_variable_names
                    except KeyboardInterrupt:
                        raise Exception("\nExiting...")

        else:
            raise Exception("Dimension (NDIM) in Optimization.cfg non existing or wrongly defined!\nExiting...")

    def objective_function_checks(self):
        """
        Function that performs checks for consistency in the definition of the objective function between the given
        configuration file and SU2 config file
        """
        # Check for consistency between Optimization.cfg and SU2.cfg
        if not self.OBJ_FUNCTION == self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION']:
            print("\nOBJECTIVE_FUNCTION in Optimization.cfg and %s" % self.SU2_CONFIG + " defined differently!")
            print("Optimization.cfg : %s" % self.OBJ_FUNCTION + "  |  %s" % self.SU2_CONFIG \
                  + " : %s\n" % self.SU2_CONFIG_IN['OBJECTIVE_FUNCTION'])
            print("OBJECTIVE_FUNCTION in %s will be changed accordingly" % self.SU2_CONFIG)
            print("Press Ctrl + C to abort")
            try:
                time.sleep(5)
                SU2_Config_change(DIR + self.SU2_CONFIG,
                                  DIR + self.SU2_CONFIG, ['OBJECTIVE_FUNCTION'], [self.OBJ_FUNCTION])
            except KeyboardInterrupt:
                raise Exception("\nExiting...")

    def restart_checks(self):
        """
        Function that stores variables and performs checks for the restarting option
        """
        # Store number of previous design iterations
        self.DesignStep_Restart = 0
        for i in range(100):
            if os.path.exists('./DESIGN/DSN_%02i/' % (i + 1)):
                self.DesignStep_Restart += 1
            elif not os.path.exists('./DESIGN/DSN_%02i/' % (i + 1)):
                break

        # Create array for relaxation factor
        # This is needed in case the restart is done where previously the RF has been changed
        # The array is updated few lines below
        self.of_relax_factor = [None] * self.DesignStep_Restart
        if self.IN['OPT_CONSTRAIN'] is not 'NONE':
            self.cons_relax_factor = [None] * self.DesignStep_Restart

        # Check if DSN directory exists
        if not os.path.exists("./DESIGN/DSN_01/"):
            raise Exception("DESIGN/DSN_** folder does not exist in the path ... could not restart optimization ...\nExiting...")

        # Check if pickle file is correctly written
        try:
            for ii in range(1, self.DesignStep_Restart + 1):
                self.OPTSET_PCKL = open("./DESIGN/DSN_%02i/optset.pickle" % ii,"rb")
                var_names = pickle.load(self.OPTSET_PCKL)                               # Variable names
                var_vals = pickle.load(self.OPTSET_PCKL)                                # Variable values
                par = pickle.load(self.OPTSET_PCKL)                                     # Design vector
                self.of_relax_factor[ii - 1] = pickle.load(self.OPTSET_PCKL)            # Optimization relaxation factor
                if self.IN['OPT_CONSTRAIN'] is not 'NONE':
                    self.cons_relax_factor[ii -1] = pickle.load(self.OPTSET_PCKL)       # Constrain relaxation factor

        except:
            raise Exception("Pickle file (*.pickle) does not exist in the DESIGN folder!\nExiting...")

        # Check consistency between pickle DVs and config. DVs
        #pdb.set_trace()
        if not var_names == self.variable_names:
            print("Variable names in Optimization.cfg do not match with variable names from previous runs")
            print("Switching RESTART to 'False' ...\n Press Ctrl + C to exit ...")
            try:
                time.sleep(5)
                self.RESTART = 'NO'
            except KeyboardInterrupt:
                raise Exception("Exiting...")

    def search_similar_mesh(self):
        """
        Function that searches amongst all the existing meshes the most similar one by calculating the RMS of the difference
        between the surface points.
        Used to perform mesh deformation on a similar mesh and avoid huge mesh deformations
        :return: DIR_MESH = Directory of the most similar mesh
        """
        # Read MoveSurface.txt for the current design step
        move_surf = np.loadtxt(self.CAD_FOLDER % self.DesignStep + self.SU2_CONFIG_IN['MOTION_FILENAME'])

        # Initialize variables
        RMS = [None] * (self.DesignStep - 1)
        # Loop through all design steps available
        step = 1
        for ii in range(self.DesignStep - 1):
            if step == 1:
                loc_comp = self.DEFORM_FOLDER % step + self.SU2_CONFIG_IN['MOTION_FILENAME']
            else:
                loc_comp = self.CAD_FOLDER % step + self.SU2_CONFIG_IN['MOTION_FILENAME']

            # Read MoveSurface.txt for the corresponding design step (variable "step")
            move_surf_comp = np.loadtxt(loc_comp)

            RMS_sum = 0
            for kk in range(len(move_surf_comp)):
                # Compute distance between points
                D = np.sqrt((move_surf[kk,1] - move_surf_comp[kk,1])**2 + (move_surf[kk,2] - move_surf_comp[kk,2])**2 + (move_surf[kk,3] - move_surf_comp[kk,3])**2)
                # Compute sum of squares for RMS
                RMS_sum += D**2
            # Compute RMS for the corresponding design step and store it in an array
            # (mesh corresponding to "step" vs. mesh corresponding to DesignStep)
            RMS[ii] = np.sqrt(RMS_sum/len(move_surf_comp))

            # Increase step
            step += 1

        # Find minimum RMS and its corresponding position
        position_min = np.argmin(RMS)

        # Locate directory with most similar mesh based on minimum RMS value
        if position_min == 0:
            DIR_MESH = DIR + self.SU2_CONFIG_IN['MESH_FILENAME']
        else:
            DIR_MESH = self.DEFORM_FOLDER % (position_min + 1) + self.SU2_CONFIG_IN['MESH_OUT_FILENAME']

        print("... Closest mesh located at: %s ..." % DIR_MESH)

        # Return directory
        return DIR_MESH

    def run_command(self, Command):
        try:
            # Executable argument added for cluster compatibility
            # If problematic, please delete it
            proc = subprocess.Popen(Command, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=subprocess.PIPE)
            return_code = proc.wait()
            message = (proc.stderr.read())
            print(message.decode('ascii'))
        except:
            print("Did subprocess crash ??? ")

    def pickleize(self,design_step_pickle):
        """
        Function that converts previous CFD and ADJ runs into pickle files (in case the run was aborted - e.g. cluster)
        :return: Pickle files of the CFD and ADJ runs
        """
        # Modules to pickleize
        pickleize_CFD = True
        pickleize_ADJ = True

        # Directories
        DIRECT_DIR  = self.DIRECT_FOLDER % design_step_pickle
        ADJ_DIR     = self.SENS_OBJ_FOLDER % (design_step_pickle, self.OBJ_FUNCTION)
        CAD_DIR     = self.CAD_FOLDER % design_step_pickle

        # Directory and file checks
        if pickleize_CFD:
            if os.path.exists(DIRECT_DIR + '/obj_fun.pickle'):
                print(
                    "obj_fun.pickle already exists in DIRECT folder and it will be removed!\nPress Ctrl + C to abort...")
                try:
                    time.sleep(5)
                    os.system('rm -rf %s' % (DIRECT_DIR + '/obj_fun.pickle'))

                    # Create new pickle file
                    os.chdir(DIRECT_DIR)
                    OBJ_PICKL = open(DIRECT_DIR + '/obj_fun.pickle', "wb")
                except KeyboardInterrupt:
                    raise Exception("Check configuration in pickleizer.py...\nExiting...")

            else:
                # Create new pickle file
                os.chdir(DIRECT_DIR)
                OBJ_PICKL = open(DIRECT_DIR + '/obj_fun.pickle', "wb")

        if pickleize_ADJ:
            if os.path.exists(ADJ_DIR + '/sens_obj.pickle'):
                print(
                    "obj_fun.pickle already exists in DIRECT folder and it will be removed!\nPress Ctrl + C to abort...")
                try:
                    time.sleep(5)
                    os.system('rm -rf %s' % (ADJ_DIR + '/sens_obj.pickle'))

                    # Create new pickle file
                    os.chdir(ADJ_DIR)
                    OBJ_SENS_PICKL = open(ADJ_DIR + '/sens_obj.pickle', "wb")
                except KeyboardInterrupt:
                    raise Exception("Check configuration in pickleizer.py...\nExiting...")

            else:
                # Create new pickle file
                os.chdir(ADJ_DIR)
                OBJ_SENS_PICKL = open(ADJ_DIR + '/sens_obj.pickle', "wb")

        else:
            raise Exception("Check options in ShapeOptimization.Pickleizer...\nExiting...")

        # Write pickle files
        if pickleize_CFD:
            # Change directory
            os.chdir(DIRECT_DIR)

            # Read history.dat
            History = np.loadtxt('history.dat', skiprows=3, delimiter=',')

            # Read converged values
            conv_values = History[-1][1:16]
            OBJ = conv_values[self.INDEX[self.OBJ_FUNCTION] - 1]

            # Write pickle file
            pickle.dump(OBJ, OBJ_PICKL)
            OBJ_PICKL.close()
            print("Created obj_fun.pickle ...")

        if pickleize_ADJ:
            # Run CAD
            os.chdir(CAD_DIR)
            self.run_cad_sens(self.CAD_CONFIG_IN)

            # Change directory
            os.chdir(ADJ_DIR)

            # Run dot_product
            Sens = self.dot_product()

            # Write sensitivities (with relax. factor) in the pickle file
            pickle.dump(Sens, OBJ_SENS_PICKL)
            OBJ_SENS_PICKL.close()
            print("Created sens_obj.pickle ...")

        print("Pickleize process ended!")

    def initialize_DVs(self):
        """
        Function to initialize design variable names
        """
        # Read variable list from configuration file
        # Thickness design varuables
        if (self.IN['DV_DICTIONARY'] == 'THICKNESS'):
            if self.CAD_CONFIG_IN['PARAMETRIZATION_TYPE'] == 'CONNECTING_ARCS':
                self.variable_names = ['dist_1', 'dist_2', 'dist_3', 'dist_4']
            else:
                self.variable_names = ['thickness_upper_1', 'thickness_upper_2', 'thickness_upper_3',
                                        'thickness_upper_4', 'thickness_upper_5', 'thickness_upper_6',
                                        'thickness_lower_1', 'thickness_lower_2', 'thickness_lower_3',
                                        'thickness_lower_4', 'thickness_lower_5', 'thickness_lower_6']

        # Section design variables
        elif (self.IN['DV_DICTIONARY'] == 'SECTION'):
            if self.CAD_CONFIG_IN['PARAMETRIZATION_TYPE'] == 'CONNECTING_ARCS':
                self.variable_names = self.blade_section_connecting_arcs
            else:
                self.variable_names = self.blade_section_camber_thickness

        # Meridional channel design variables
        elif (self.IN['DV_DICTIONARY'] == 'CHANNEL'):
            self.variable_names = self.meridional_channel_names

        # All design variables
        elif (self.IN['DV_DICTIONARY'] == 'ALL'):
            if self.CAD_CONFIG_IN['PARAMETRIZATION_TYPE'] == 'CONNECTING_ARCS':
                self.variable_names = self.blade_section_connecting_arcs + self.meridional_channel_names
            else:
                self.variable_names = self.blade_section_camber_thickness + self.meridional_channel_names
        else:
            self.variable_names = self.IN['DV_DICTIONARY']

        if isinstance(self.variable_names, str):
            self.variable_names = [self.variable_names]                     # Store single DV in a list

        # self.design_variable_checks()                                     # Perform consistency checks on DVs

        self.N_variables = len(self.variable_names)                         # Total number of design variables defined in config

        self.N_DVs = 0                                                      # Number of DVs (counting for variable names
        for name in self.variable_names:                                    # that might have more than one value)
            self.N_DVs += len(self.CAD_CONFIG_IN[name])

    def runtime_options(self):
        """
        Initializes and performs checks on several options that can be defined in the configuration file
        """

        if self.OPERATION_TYPE != 'PLOT_VALIDATION':
                self.SU2_CONFIG_IN  = ReadUserInput(DIR + self.SU2_CONFIG)      # Read SU2 configuration file
                self.CAD_CONFIG_IN  = ReadUserInput(DIR + self.CAD_CONFIG)      # Read CAD configuration file

                self.OBJ_FUNCTION = self.IN['OBJECTIVE_FUNCTION']               # Objective function defined by the user
                self.objective_function_checks()                                # Perform consistency checks on OBJ fun
        else:
            self.CAD_CONFIG_IN = ReadUserInput(DIR + self.CAD_CONFIG)

        self.initialize_DVs()
        #try:
        self.RESTART = self.IN['RESTART']                                   # Restart Optimization from existing DSN_** folders
        if self.RESTART == 'YES':
            self.restart_checks()                                           # Perform checks for restarting
        #except:
         #   self.RESTART = 'NO'

        try:
            self.SKIP_REF_CFD = self.IN['SKIP_REF_CFD']                         # Skip reference CFD for finite-differences (requires existing disc. adjoint run)
        except:
            self.SKIP_REF_CFD = 'NO'

        try:
            self.SKIP_ADJOINT = self.IN['SKIP_ADJOINT']                         # Skip adjoint calculation for disc. adjoint and runs only dot_product
        except:
            self.SKIP_ADJOINT = 'NO'

        try:
            self.REMESH = self.IN['REMESH']                                     # Remesh
            if self.REMESH == 'YES':
                self.REMESH = True
                self.MESH_SEARCH = False
            else:
                self.REMESH = False
        except:
            self.REMESH = False

        try:
            self.REMESH_STEP = int(self.IN['REMESH_STEP'][0])                   # Remesh step as defined in config. file
        except:
            pass

BladeOpti = ShapeOptimization(IN)
BladeOpti.start_process()
