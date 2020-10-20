###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

############################### FILE NAME: blade_output.py ####################################
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

#----------------------------------------------------------------------------------------------------------------------#
# Define BladeOutput class
#----------------------------------------------------------------------------------------------------------------------#
class BladeOutput:
    def __init__(self, blade_in):

        # Declare blade_in variables as instance variables
        self.blade_in = blade_in

        # Configuration options
        self.OPERATION_TYPE = blade_in.OPERATION_TYPE
        self.NDIM = blade_in.NDIM
        self.N_SECTIONS = blade_in.N_SECTIONS
        self.N_points = blade_in.N_points

        # Blade geometry
        self.u = blade_in.u
        self.v = blade_in.v
        self.surface_coordinates = np.real(blade_in.surface_coordinates)
        self.surface_normals = np.real(blade_in.surface_normals)
        self.surface_sensitivity = blade_in.surface_sensitivity

        # Meridional channel geometry
        self.u_hub = blade_in.u_hub
        self.u_shroud = blade_in.u_shroud
        self.hub_coordinates = np.real(blade_in.hub_coordinates)
        self.shroud_coordinates = np.real(blade_in.shroud_coordinates)

        # Create output directory
        self.full_path = None
        os.system("rm -rf output")
        if os.path.exists("output") is False:
            os.mkdir("output")


    # ---------------------------------------------------------------------------------------------------------------- #
    # Print output files with the blade geometry
    # ---------------------------------------------------------------------------------------------------------------- #
    def print_blade_coordinates(self, filename='surface_coordinates', path=os.getcwd()):

        """ Print surface coordinates in a .csv file """

        print("Writting blade surface coordinates...", end='                  ')

        # Create output directory in case it does not exist
        full_path = path + '/output/coordinates/'
        if os.path.exists(full_path) is False:
            os.mkdir(full_path)

        file = open(full_path + filename + '.csv', 'w')

        if self.NDIM == 3:
            header = "\"Point\",\t\"x_coord\",\t\"y_coord\",\t\"z_coord\",\t\"u\",\t\"v\"\n"
            file.write(header)
            for i in range(self.N_points):
                file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%.16f,\t%.16f\n' %
                           (i,                                              # Index
                            self.surface_coordinates[2, i],                 # SU2 x coordinate (radial)
                            self.surface_coordinates[1, i],                 # SU2 y coordinate (tangential)
                            self.surface_coordinates[0, i],                 # SU2 z coordinate (axial)
                            self.u[i], self.v[i]))                          # Index and (u,v) pair
        elif self.NDIM == 2:
            header = "\"Point\",\t\"x_coord\",\t\"y_coord\",\t\"u\",\t\"v\"\n"
            file.write(header)
            for i in range(self.N_points):
                file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%.16f\n' %
                           (i,                                              # Index
                            self.surface_coordinates[0, i],                 # x coordinate (axial)
                            self.surface_coordinates[1, i],                 # y coordinate (tangential)
                            self.u[i], self.v[i]))                          # Index and (u,v) pair
        else:
            raise Exception("The number of dimensions must be 2 or 3")

        file.close()

        print("Done!")


    def print_hub_coordinates(self, filename='hub_coordinates', path=os.getcwd()):

        """ Print hub coordinates in a .csv file """

        print("Writing hub surface coordinates...", end='                    ')

        # Create output directory in case it does not exist
        full_path = path + '/output/coordinates/'
        if os.path.exists(full_path) is False:
            os.mkdir(full_path)

        # Print the hub surface coordinates
        file = open(full_path + filename + '.csv', 'w')
        header = "\"Point\",\t\"x_coord\",\t\"y_coord\",\t\"z_coord\",\t\"u\"\n"
        file.write(header)
        for i in range(np.shape(self.hub_coordinates)[1]):
            file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%.16f\n'
                       % (i,                                            # Index
                          self.hub_coordinates[1, i],                   # SU2 x coordinate (radial)
                          0.00,                                         # SU2 y coordinate (tangential)
                          self.hub_coordinates[0, i],                   # SU2 z coordinate (axial)
                          self.u_hub[i]))                               # u-parameter
        file.close()
        print("Done!")


    def print_shroud_coordinates(self, filename='shroud_coordinates', path=os.getcwd()):

        """ Print shroud coordinates in a .csv file """

        print("Writting shroud surface coordinates...", end='                 ')

        # Create output directory in case it does not exist
        full_path = path + '/output/coordinates/'
        if os.path.exists(full_path) is False:
            os.mkdir(full_path)

        # Print the shroud surface coordinates
        file = open(full_path + filename + '.csv', 'w')
        header = "\"Point\",\t\"x_coord\",\t\"y_coord\",\t\"z_coord\",\t\"u\"\n"
        file.write(header)
        for i in range(np.shape(self.shroud_coordinates)[1]):
            file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%.16f\n'
                       % (i,                                            # Index
                          self.shroud_coordinates[1, i],                # SU2 x coordinate (radial)
                          0.00,                                         # SU2 y coordinate (tangential)
                          self.shroud_coordinates[0, i],                # SU2 z coordinate (axial)
                          self.u_shroud[i]))                            # u-parameter
        file.close()
        print("Done!")


    def print_sensitivity(self, filename='grad', path=os.getcwd()):

        """ Print surface sensitivities in .csv files """

        print("Writting surface sensitivity...", end='                        ')

        if self.OPERATION_TYPE == 'SENSITIVITY':

            # Create output directory in case it does not exist
            full_path = path + '/output/sensitivities/'
            if os.path.exists(full_path) is False:
                os.mkdir(full_path)

            for key in self.surface_sensitivity.keys():
                file = open(full_path + filename + '_' + key + '.csv', 'w')
                if self.NDIM == 3:
                    header = "\"Point\"\t\"x_coord\"\t\"y_coord\"\t\"z_coord\"\t\" x_sens\"\t\"y_sens\"\t\"z_sens\"\t\" x_normal\"\t\"y_normal\"\t\"z_normal\"\n"
                    file.write(header)
                    for i in range(self.N_points):
                        file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e\n'
                                   % (i,                                        # Index
                                      self.surface_coordinates[2, i],           # SU2 x coordinate (radial)
                                      self.surface_coordinates[1, i],           # SU2 y coordinate (tangential)
                                      self.surface_coordinates[0, i],           # SU2 z coordinate (axial)
                                      self.surface_sensitivity[key][2, i],      # SU2 x sensitivity
                                      self.surface_sensitivity[key][1, i],      # SU2 y sensitivity
                                      self.surface_sensitivity[key][0, i],      # SU2 z sensitivity
                                      self.surface_normals[2, i],               # SU2 x normal
                                      self.surface_normals[1, i],               # SU2 y normal
                                      self.surface_normals[0, i]))              # SU2 z normal

                elif self.NDIM == 2:
                    header = "\"Point\"\t\"x_coord\"\t\"y_coord\"\t\" x_sens\"\t\"y_sens\"\t\" x_normal\"\t\"y_normal\"\n"
                    file.write(header)
                    for i in range(self.N_points):
                        file.write('%i,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e,\t%+.16e\n'
                                   % (i,                                        # Index
                                      self.surface_coordinates[0, i],           # x coordinate (axial)
                                      self.surface_coordinates[1, i],           # y coordinate (tangential)
                                      self.surface_sensitivity[key][0, i],      # x sensitivity
                                      self.surface_sensitivity[key][1, i],      # y sensitivity
                                      self.surface_normals[0, i],               # x normal
                                      self.surface_normals[1, i]))              # y normal
                else:
                    raise Exception("The number of dimensions must be 2 or 3")

                file.close()

        print("Done!")


    # ---------------------------------------------------------------------------------------------------------------- #
    # Write mesh generation files
    # ---------------------------------------------------------------------------------------------------------------- #
    def write_mesh_file(self, path=os.getcwd()):

        """ Write a mesh file that can be used to generate a CFD grid

            - Meshes in 2D are generated with UMG2
            - Meshes in 3D are generated with ANSYS Turbogrid

        """

        # Create output directory in case it does not exist
        self.full_path = path + '/output/mesh_files/'
        if os.path.exists(self.full_path) is False:
            os.mkdir(self.full_path)

        # Write the 2D or 3D mesh generation files
        if self.NDIM == 2:
            self.make_UMG2_meshfile()
        elif self.NDIM == 3:
            self.make_TurboGrid_meshfile()
        else:
            raise Exception('The number of dimensions must be "2" or "3"')


    def make_UMG2_meshfile(self):

        """ Write file to generate an mesh using UMG2 """

        print("Writing UMG2 mesh generation files...", end='                  ')

        blade = open(self.full_path + "/blade.geometry", 'w')
        blade.write("Not yet implemented")  # TODO implement 2D mesh generation

        print('Done!')


    def make_TurboGrid_meshfile(self):

        """ Write file to generate an mesh using ANSYS TurboGrid  """

        print("Writing TurboGrid mesh generation files...", end='               ')

        # Define the number of blade sections and the number of points per section
        n_sections = self.N_SECTIONS
        n_points = 500

        # Define the (u,v) parametrization (with an offset ´h´ from limits)
        h = 1e-5
        v = np.linspace(0.00 + h, 1.00 - h, n_sections)
        u = np.linspace(0.00 + h, 1.00 - h, n_points)

        # Write blade.crv file
        blade_file = open(self.full_path + "blade.crv", 'w')
        section_header = "#blade_%i\n"
        index = 1
        for v_section in v:
            blade_file.write(section_header % index)
            section_coordinates = np.real(self.blade_in.get_surface_coordinates(u=u, v=np.ones(n_points)*v_section))
            for i in range(n_points):
                blade_file.write("%+.16e\t%+.16e\t%+.16e\n" %
                                 (section_coordinates[2, i],            # SU2 x coordinate (radial)
                                  section_coordinates[1, i],            # SU2 y coordinate (tangential)
                                  section_coordinates[0, i]))           # SU2 z coordinate (axial)
            index = index+1
        blade_file.close()

        # Write hub.crv file
        hub_file = open(self.full_path + "hub.crv",'w')
        hub_coordinates = np.real(self.blade_in.get_extended_hub_coordinates(u))
        for i in range(n_points):
            hub_file.write("%+.16e\t%+.16e\t%+.16e\n" % (hub_coordinates[1, i], 0.00, hub_coordinates[0, i]))
        hub_file.close()

        # Write shroud.crv file
        shroud_file = open(self.full_path + "shroud.crv",'w')
        shroud_coordinates = np.real(self.blade_in.get_extended_shroud_coordinates(u))
        for i in range(n_points):
            shroud_file.write("%+.16e\t%+.16e\t%+.16e\n" % (shroud_coordinates[1, i], 0.00, shroud_coordinates[0, i]))
        shroud_file.close()

        print('Done!')
