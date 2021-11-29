#!/usr/bin/env python3
###############################################################################################
#                    ____                 ____  _           _                                 #
#                   |  _ \ __ _ _ __ __ _| __ )| | __ _  __| | ___                            #
#                   | |_) / _` | '__/ _` |  _ \| |/ _` |/ _` |/ _ \                           #
#                   |  __/ (_| | | | (_| | |_) | | (_| | (_| |  __/                           #
#                   |_|   \__,_|_|  \__,_|____/|_|\__,_|\__,_|\___|                           #
#                                                                                             #
###############################################################################################

################################# FILE NAME: PlotBlade.py #####################################
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
import numpy as np
import matplotlib.pyplot as plt
import pdb
import sys
import os
import time

#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
sys.path.append(os.getcwd() + '/../parablade')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from parablade.common import PrintBanner
from parablade.common.config import *
from parablade.blade_3D import Blade3D
from parablade.blade_plot import BladePlot
from parablade.blade_output import BladeOutput

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
print('The blade geometry was generated in %(my_time).5f seconds' % {'my_time': time.time() - t})

#---------------------------------------------------------------------------------------------#
# Plot the blade profile
#---------------------------------------------------------------------------------------------#
plot_blade_object = BladePlot(blade_object)
plot_blade_object.make_plots()

#---------------------------------------------------------------------------------------------#
# End of file
#---------------------------------------------------------------------------------------------#