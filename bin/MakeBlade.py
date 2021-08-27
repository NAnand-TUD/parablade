#!/usr/bin/env python3
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
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pdb
import sys
import os
import time

#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common.common import *
from common.config import *
from src.blade_3D import Blade3D
from src.blade_plot import BladePlot
from src.blade_output import BladeOutput

#---------------------------------------------------------------------------------------------#
# Print ParaBlade Banner
#---------------------------------------------------------------------------------------------#
PrintBanner()

#---------------------------------------------------------------------------------------------#
# Initializing the DIR and load the configuration file
#---------------------------------------------------------------------------------------------#
DIR = os.getcwd() + '/'

try:
    INFile = DIR + sys.argv[-1]
except:
    INFile = DIR + 'blade.cfg'      # Default File name

try:
    IN = ReadUserInput(INFile)
except:
    raise Exception('\n\n\n''Something went wrong when reading the configuration file,exiting the program...'
                    '\n\nTo call MakeBlade.py from terminal type:'
                    '\n\tMakeBlade.py <configuration file name>')

#---------------------------------------------------------------------------------------------#
# Generation of blade profile
#---------------------------------------------------------------------------------------------#
t = time.time()
blade_object = Blade3D(IN)
blade_object.make_blade()
print('The blade geometry was generated in %(my_time).5f seconds\n' % {'my_time': time.time() - t})

# --------------------------------------------------------------------------------------------#
# Write the output files
# --------------------------------------------------------------------------------------------#
blade_output_object = BladeOutput(blade_object)
blade_output_object.print_blade_coordinates()
blade_output_object.print_hub_coordinates()
blade_output_object.print_shroud_coordinates()
if IN['OPERATION_TYPE'] == "SENSITIVITY":
    blade_output_object.print_sensitivity()

# --------------------------------------------------------------------------------------------#
# Write Mesh files
# --------------------------------------------------------------------------------------------#
blade_output_object.write_mesh_file()

#---------------------------------------------------------------------------------------------#
# End of file
#---------------------------------------------------------------------------------------------#