#!/usr/bin/env python3

""" Regression test for the thickness computation

    Compute thickness measures for 2D NACA 4 digit series airfoils analytically, numerically with the
    thickness calculation methods and the numerical thickness of the matched airfoil.
    
    With the latter, once known the design parameters, also calculate their sensitivity 
    using different differentiation methods:

        1) Forward finite differences
        2) Central finite differences
        3) Complex step

    Author: Ricardo Puente
    Date: 10/2019

"""

import os
import sys
import numpy as np
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm 
import inspect

#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common import PrintBanner
from config import *
from blade_3D import Blade3D
from blade_match import BladeMatch
from blade_output import BladeOutput
from CAD_functions import split_curve

# Rebuild the matplotlib font cache
fm._rebuild()
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['savefig.dpi'] = 600
font = {'weight' : 'normal',
        'size'   : 20}

plt.rc('font', **font)

def NACA_4_digit(t,m,p,c,Npoints=500):
    
    """ 

        Parameters
        ----------
        Npoints: Number of axial points
        t      : Maximum thickness as fraction of the chord
        m      : Maximum camber as fraction of the chord
        p      : Location of maximum camber
        c      : Chord
            
        Returns
        ----------
        x_up : x upper side coordinate array
        y_up : y upper side coordinate array
        x_dw : x lower side coordinate array
        y_dw : y lower side coordinate array
        xc   : x camber coordinate array
        yc   : y lower side coordinate array
        th   : Thickness distribution
        dth  : Thickness distribution derivative
        x_mxth : Location of the maximum thickness
                    
        Author: Ricardo Puente, 0/2020
                r.puente@imperial.ac.uk

    """
    
    # Basic input sanity checks
    if p < 0.:
        p = 0.0001
    elif p > 1.:
        p = 1.
    elif p == 0.:
        raise ValueError('Maximum camber location p cannot be 0')
        
    if c < 0.:
        c = 1.
    
    if m<0.:
        m = 0.
    
    x = np.linspace(1.e-8,1.,Npoints)
    x_up = np.zeros(Npoints)
    x_dw = np.zeros(Npoints)
    y_up = np.zeros(Npoints)
    y_dw = np.zeros(Npoints)
    yc   = np.zeros(Npoints)
    
    
    # Symmetric part (thickness)
    alpha = 0.2969
    beta = -0.126
    gamma = -0.3516
    delta = 0.2843
    kappa = -0.1015
    
    x2 = np.power(x,2.)
    x3 = np.power(x,3.)
    x4 = np.power(x,4.)
    
    th_half  = 5.*t*(alpha*np.sqrt(x) + beta*x + gamma*x2 + delta*x3 + kappa*x4)
    dth_half = 5.*t*(0.5*alpha/np.sqrt(x) + beta + 2.*gamma*x + 3.*delta*x2 + 4.*kappa*x3)
    
    m_p2  = m/(p*p)
    m_1p2 = m/((1.-p)*(1.-p))
    for i,xi in enumerate(x):
        if xi < p:
            yc[i]  = m_p2*(2.*p*xi-xi*xi)
            dyc = 2.*m_p2*(p-xi)
        else:
            yc[i] =m_1p2*((1.-2.*p)+2.*p*xi-xi*xi)
            dyc = 2.*m_1p2*(p-xi)
            
        theta = np.arctan(dyc)
        x_up[i] = xi-th_half[i]*np.sin(theta)
        x_dw[i] = xi+th_half[i]*np.sin(theta)
        
        y_up[i] = yc[i]+th_half[i]*np.cos(theta)
        y_dw[i] = yc[i]-th_half[i]*np.cos(theta)
        
    
    # Scale with the chord
    th=2.*th_half*c
    dth=2.*dth_half*c
    x_up=x_up*c
    x_dw=x_dw*c
    y_up=y_up*c
    y_dw=y_dw*c
    
    # This parameter is a fixed value in the NACA 4 digit series
    x_mxth = c*0.3
    
    xc = x*c
    yc = yc*c
        
    
    return x_up,y_up,x_dw,y_dw,xc,yc,th,dth,x_mxth


def AirfoilMatch(INFile,airfoilToMatch,filesExist=False):

    # Write the airforil to be matched in a file    
    coordsFileName = "naca.txt"
    
    coordsFile = open(coordsFileName,'w')
    for i in range(len(airfoilToMatch[0])):
        line = str(i+1)+"\t"+str(airfoilToMatch[0][i])+"\t"+str(airfoilToMatch[1][i])+"\n"
        coordsFile.write(line)

    coordsFile.close()
    
    # Create BladeMatch object
    if filesExist==False:
        options = {'view_xy'            : 'no',    # 2D Recommended
                   'view_xR'            : 'no',     # 3D Recommended
                   'view_yz'            : 'no',     # 3D Optional
                   'view_3D'            : 'no',     # 3D Recommended
                   'error_distribution' : 'no'}
        
        IN = ReadUserInput(INFile)
        matched_blade_object = BladeMatch(IN, coarseness=1, plot_options=options)

        # Match the blade automatically (optimization mode)
        matched_blade_object.match_blade(matching_mode='DVs')
    
    # Write the output config file
    OUT =  "matched_parametrization.cfg"
    path=os.getcwd()
    full_path = path + '/output_matching/'
    ifi =  full_path + OUT 
    
    # Create a blade obect with the matched parameters to return
    try:
       INnew = ReadUserInput(ifi)
    except:
        raise Exception('\n\n\n''Problems in AirfoilMatch(): Something went wrong when reading the matched airfoil configuration file!')
        
    matched_airfoil = Blade3D(INnew)
    
    return matched_airfoil
    

def NACA_validation(t,m,p,c,file,save_plots=False,Npoints=500):
    
    #########################
    # Generate a NACA profile
    #########################


    x_up,y_up,x_dw,y_dw,xc,yc,th,dth,x_mxth = NACA_4_digit(t,m,p,c,Npoints)

    #####################################
    # Compute the camber line numerically
    ######################################
    
    # Check for strictly increasing sets of points or the thickness calculation will fail
    def checkIfStriclyIncreasing(L):
        return all(x<y for x, y in zip(L, L[1:]))
    
    def fixCoordinates(x,y):
    
        xs=x.copy()
        ys=y.copy()
        if checkIfStriclyIncreasing(x)==False:
            xs = x[1:]
            ys = y[1:]
            fixCoordinates(xs,ys)
            
        return xs,ys
        
       
   
    x_up,y_up = fixCoordinates(x_up,y_up)
    x_dw,y_dw = fixCoordinates(x_dw,y_dw)  
    
    
    section_upper_side_list = []
    section_upper_side_list.append(x_up)
    section_upper_side_list.append(y_up)
    section_lower_side_list = []
    section_lower_side_list.append(x_dw)
    section_lower_side_list.append(y_dw)

    section_upper_side = np.array(section_upper_side_list)
    section_lower_side = np.array(section_lower_side_list)
    section_thickness_dist,max_th,lack_of_monotonicity,camber,section_thickness_dist_der = Blade3D.get_section_thickness_properties(section_upper_side,section_lower_side)  
 
    #####################################
    # Match airfoil shape
    ######################################
    
    # Build an array with the airfoil coordinates
    INfile = "Template.cfg"
    xtotal = np.concatenate((x_up,np.flip(x_dw)), axis=0)
    ytotal = np.concatenate((y_up,np.flip(y_dw)), axis=0)
    naca_coords = np.zeros((2,len(xtotal)),dtype=float)
    naca_coords[0]=xtotal
    naca_coords[1]=ytotal
    airfoilMatch = AirfoilMatch(INfile,naca_coords)
    match_coords = airfoilMatch.section_coordinates[0]
    upper_side_M,lower_side_M = split_curve(match_coords)
    section_thickness_dist_M,max_th_M,lack_of_monotonicity_M,camber_M,section_thickness_dist_der_M = airfoilMatch.get_section_thickness_properties(upper_side_M,lower_side_M)
    
    ##############################################
    # Print the maximum thickness and its location
    ###############################################
    max_th_index = np.where(section_thickness_dist[1]==max_th)

    out = str([t,m,p]) + "\t" + str(t) + "\t" + str(max_th.real) + "\t" + str(x_mxth) + "\t" \
        + str(section_thickness_dist[0][max_th_index].real) + "\t" + str(lack_of_monotonicity.real) + "\n" 
    file.write(out)

    ######################
    # Compute errors
    #######################

    # First interpolate the analytical thickness in the basis of the numerical one
    InterpolateAnalytic = interpolate.CubicSpline(xc,th)
    errorDirect = np.zeros(len(section_thickness_dist[0]))
    
    for i,th_i in enumerate(section_thickness_dist[1]):
        xi = section_thickness_dist[0][i]
        evalTh = InterpolateAnalytic(xi)
        errorDirect[i] = abs(th_i-evalTh)/evalTh
        
    errorDirect_M = np.zeros(len(section_thickness_dist_M[0]))
    
    for i,th_i in enumerate(section_thickness_dist_M[1]):
        xi = section_thickness_dist_M[0][i]
        evalTh = InterpolateAnalytic(xi)
        errorDirect_M[i] = abs(th_i-evalTh)/evalTh


    ##################
    # Generate plots 
    ##################


    fig, (ax1, ax2, ax4) = plt.subplots(3,figsize=(10,10),sharex=True)
    ax1.set_ylabel('$y$', color='k', labelpad=12)
    ax4.set_xlabel('$x$', color='k', labelpad=12)
    ax2.set_ylabel('$Th$', color='r', labelpad=12)
    ax4.set_ylabel(r'$\frac{|Th(x)-Th_{Num}(x)|}{Th(x)}$', color='k', labelpad=12)

    #  Airfoil plots
    line, = ax1.plot(x_up, y_up)
    line.set_linewidth(1.25)
    line.set_linestyle("None")
    line.set_marker(".")
    line.set_color("k")
            
    line, = ax1.plot(x_dw, y_dw)
    line.set_linewidth(1.25)
    line.set_linestyle("None")
    line.set_marker(".")
    line.set_color("k")
    line.set_label('Airfoil')
        
    line, = ax1.plot(xc, yc)
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    line.set_color("b")
    line.set_label('Analytical camber')

    line, = ax1.plot(camber[0],camber[1])
    line.set_linewidth(1.25)
    line.set_linestyle("None")
    line.set_marker("P")
    line.set_color("r")
    line.set_label('Numerical camber')

    line, = ax1.plot(camber_M[0],camber_M[1])
    line.set_linewidth(1.25)
    line.set_linestyle("None")
    line.set_marker("X")
    line.set_color("g")
    line.set_label('Matched camber')
    
    line, = ax1.plot(match_coords[0],match_coords[1])
    line.set_linewidth(1.25)
    line.set_linestyle("-.")
    line.set_color("k")
    line.set_label('Matched airfoil')

    # Thickness and thickness derivative plot
    line, = ax2.plot(xc,th)
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    line.set_color("r")
    line.set_label('Analytical')


    ax3 = ax2.twinx()
    ax3.set_ylabel(r'$\frac{d Th}{d x}$', color='b', labelpad=12)
    ax3.set_ylim([min(dth),5])
    line,= ax3.plot(xc,dth)
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    line.set_color("b")
    line.set_label('Analytical')


    line, = ax2.plot(section_thickness_dist[0],section_thickness_dist[1])
    line.set_linewidth(1.25)
    line.set_linestyle("-.")
    line.set_color("r")
    line.set_label('Numerical')


    line,= ax3.plot(section_thickness_dist[0],section_thickness_dist_der)
    line.set_linewidth(1.25)
    line.set_linestyle("-.")
    line.set_color("b")
    line.set_label('Numerical')


    line, = ax2.plot(section_thickness_dist_M[0],section_thickness_dist_M[1])
    line.set_linewidth(1.25)
    line.set_linestyle(":")
    line.set_color("r")
    line.set_label('Matched')
    
    line,= ax3.plot(section_thickness_dist_M[0],section_thickness_dist_der_M)
    line.set_linewidth(1.25)
    line.set_linestyle(":")
    line.set_color("b")
    line.set_label('Matched')
    
    
    # Thickness error plot
    ax4.set_ylim([0.,0.5])
    line, = ax4.plot(camber[0],errorDirect)
    line.set_linewidth(1.25)
    line.set_linestyle("-")
    line.set_color("r")
    line.set_label('Numerical')
    
    line, = ax4.plot(camber_M[0],errorDirect_M)
    line.set_linewidth(1.25)
    line.set_linestyle(":")
    line.set_color("r")
    line.set_label('Matched')


    ax1.legend(frameon=False,bbox_to_anchor=(1.,1.))
    ax2.legend(frameon=False,bbox_to_anchor=(-0.15,0.75))
    ax3.legend(frameon=False,bbox_to_anchor=(1.5,0.75))
    ax4.legend(frameon=False,bbox_to_anchor=(1.5,0.75))

        
    if save_plots:
        plt.savefig(figname,format=self.figFormat,dpi=self.dpi)
    else:
        plt.show()
    

def main():
    #---------------------------------------------------------------------------------------------#
    # Print ParaBlade Banner
    #---------------------------------------------------------------------------------------------#
    PrintBanner()
    
    #---------------------------------------------------------------------------------------------#
    # Naca airfoil parameters
    #---------------------------------------------------------------------------------------------#
    t = [ 0.05, 0.1,0.2]
    m = [ 0.0,0.1,0.2]
    p = [0.2,0.5,0.8]
    # t = [0.1]
    # m = [0.0]
    # p = [0.2]
    c = 1.
    
    save_plots = False
    
    file = open("thickness_regression_output.txt","w")
    header = "Case (t,m,p)\tMax Th Analytical \t Max Th Numerical \t X(max_th) Analytical \t X(max_th) Numerical \t Lack of monotonicity \n"
    file.write(header)
    
    try:
        for ti in t:
            for mi in m:
                for pi in p:
                    NACA_validation(ti,mi,pi,c,file,save_plots)
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        message = "Raised in %s" % inspect.trace()
        print(message)
                
    file.close()


if __name__== "__main__":
  main()