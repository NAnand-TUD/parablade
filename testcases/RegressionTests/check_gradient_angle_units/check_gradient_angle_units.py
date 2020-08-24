#!/usr/bin/env python3

# Importing general packages
import numpy as np
import matplotlib.pyplot as plt
import pdb
import sys
import os

# Setting Environment
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

# Importing ParaBlade classes and functions
from common import *
from config import *
from blade_3D import Blade3D
from blade_plot import BladePlot
from blade_output import BladeOutput

# Define (u,v) parametrization
u = np.linspace(0.00, 1.00, 10)
v = 0.50 + 0.00 * u

# Compute the coodinates and sensitivity of the base blade
IN = ReadUserInput('blade_1.cfg')
blade_object = Blade3D(IN)
blade_object.make_blade()
C1 = blade_object.get_surface_coordinates(u, v)
S1 = blade_object.get_surface_sensitivity(u, v)["stagger_0"]

# Compute the coordinates of the perturbed blade
IN = ReadUserInput('blade_2.cfg')
blade_object = Blade3D(IN)
blade_object.make_blade()
C2 = blade_object.get_surface_coordinates(u, v)

# Compute the sensitivity by finite differences on the base and perturbed blades
S2 = (C2 - C1) / 0.01

# Print the results
for (x1, y1, z1), (x2, y2, z2)  in zip(np.real(S1.transpose()), np.real(S2.transpose())):
    print('%+.6e %+.6e %+.6e \t %+.6e %+.6e %+.6e \t %+.6e %+.6e %+.6e' % (x1, x2, x2-x1, y1, y2, y2-y1, z1, z2, z2-z1))