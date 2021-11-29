#!/usr/bin/python3
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
# author: Pablo                                                                               |
#    :MSc Candidates,                                                                         |
#    :Power and Propulsion                                                                    |
#    :TU Delft,                                                                               |
#    :The Netherlands                                                                         |
#                                                                                             |
#                                                                                             |
# Description:                                                                                |
#                                                                                             |
#=============================================================================================#

from numpy import * 
from numpy.linalg import inv
from matplotlib import pyplot as plt
from cycler import cycler
import os
import numpy as np
import math

#import plt_config
import optparse
#plt.rcParams["font.family"]="Times_New_Roman_Normal.ttf"
#csfont={'fontname':'Times New Roman'}
import pdb

parser = optparse.OptionParser()
parser.add_option('-a', help='Adjoint file location', dest='AD')
parser.add_option('-f', help='Finite difference file location', dest='FD')
parser.add_option('-o', help='Objective func ENTG, KE, PL, CAD', dest='header')
parser.add_option('-n', help='Name of the eps file to be written', dest='name')
parser.add_option('-s', help='Show the image yes to show no to save', dest='save',default='no')
(opts, args) = parser.parse_args()
##############################################################################################
#   INPUTS                                                                                   #
##############################################################################################
# Input file converted from Tecplot format to CSV
Location = os.getcwd()

fac = 1

AD = Location + '/' + opts.AD
FD = Location + '/' + opts.FD

if opts.header == 'ENTG':
        index = 3
elif opts.header == 'KE':
        index =2
elif opts.header == 'PL':
        index =1
elif opts.header == 'FD':
        index =-2
        fac = -1
else:
        print ("ENTG or KE")
        
        
AD_Loaded = np.genfromtxt(AD,dtype='float',skip_header=1,delimiter = ',')
FD_Loaded = np.genfromtxt(FD,dtype='float',skip_header=1,delimiter = ',')

plt.plot(AD_Loaded[:,1*fac],FD_Loaded[:,index],'k*')
#plt.plot(AD_Loaded[:,1],FD_Loaded[:,index],'ko')
[x_min,x_max,y_min,y_max] = plt.axis('equal')
x,y = np.meshgrid(np.linspace(x_min*1.1,x_max*1.1,101),np.linspace(y_min*1.1,y_max*1.1,101))
axes = plt.gca()
[y_max,y_min] = axes.get_ylim()
z = abs(np.arctan(1)-np.arctan(abs(y/x)))*100
cnt = axes.contourf(x,y,z,50)
plt.colorbar(cnt)

#plt.plot([y_min,y_max],[y_min,y_max],'g--',label=r'AD == FD',linewidth = 1.5)
#plt.plot([y_min,y_max],[y_min*0.95,y_max*0.95],'b--',label=r'5% Error',linewidth = 1.0)
#plt.plot([y_min*0.95,y_max*0.95],[y_min,y_max],'b--',linewidth = 1.0)
#plt.plot([y_min,y_max],[y_min*0.9,y_max*0.9],'y--',label=r'10% Error',linewidth = 0.5)
#plt.plot([y_min*0.9,y_max*0.9],[y_min,y_max],'y--',linewidth = 0.5)
#plt.plot([y_min,y_max],[y_min*0.8,y_max*0.8],'r--',label=r'20% Error',linewidth = 0.25)
#plt.plot([y_min*0.8,y_max*0.8],[y_min,y_max],'r--',linewidth = 0.25)
plt.plot(AD_Loaded[:,1*fac],FD_Loaded[:,index],'r*',label=r'DVs')
#plt.plot([y_min*0.9,y_max*0.9],[y_min*0.9,y_max*0.9],'y--',label=r'AD != FD')
plt.title(r"FD vs AD Gradient Validation")
plt.xlabel(r"$\frac{dJ}{d\alpha}_{AD}$")
plt.ylabel(r"$\frac{dJ}{d\alpha}_{FD}$")
plt.legend()
plt.axis('equal')
plt.grid('off')
pdb.set_trace()
if opts.save == 'no':
        plt.show()
else:
        plt.savefig(opts.name+'.eps',bbox_inches="tight",format='eps', dpi=2000)
