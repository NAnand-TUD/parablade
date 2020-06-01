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

import numpy as np
import re
import pdb
import copy
import cmath

# TODO change name convention, functions should be all lower case with underscores for readibility if necesary
# TODO no abbreviations!
def WriteConfigFile(OUTFile,IN):
    #pdb.set_trace()
    for key in IN:
        input = str(IN[key])
        output = input.replace('[','')
        output1 = output.replace(']', '')
        try:
            OUTFile.write("%s=%f\n"%(key,np.real(float(output1))))
        except:
            #pdb.set_trace()
            OUTFile.write("%s=%s\n" % (key, output1))

def SU2_Config_change(name,dest,par,vals):
    '''
    Changes parameter for SU2 configuration file.
    Usage SU2_Config_change(<infile abs dest>,<outfile abs dest>,[<parameter to be changed seperated by ,>],[<values separated by commas])>])
    '''
    IN = {}
    infile = open(name,'r')
    for line in infile:
      words = re.split(r'=|%|\n|#',line)
      if not any(words[0] in s for s in ['\n','%',' ','#']):
        words = list(filter(None,words))
        IN[words[0]] = words[1]
    #changing values
    for i in range(len(par)):
        IN[par[i]]=str(vals[i])
    outfile=open(dest,'w')
    for key, value in IN.items():
        if key not in ('DV_VALUE' , 'DV_PARAM'):
            outfile.write("%s= %s\n"%(key,value))
    outfile.write("%s= %s\n" % ('DV_PARAM', IN['DV_PARAM']))
    outfile.write("%s= %s\n" % ('DV_VALUE', IN['DV_VALUE']))
    #outfile.write('DV_KIND= HICKS_HENNE\nDV_PARAM= ( 0.0, 0.05)\nDEFINITION_DV= ( 1 , 1.0 | wall1,wall2  | 0.0 , 0.05  )\nNUMBER_PART= 8\nWRT_CSV_SOL= YES\nGRADIENT_METHOD= DISCRETE_ADJOINT\nDV_MARKER= ( WING )\nDV_VALUE= 0.001')

def WriteSU2ConfigFile(ConfigName,TYPE):
    with open(ConfigName, 'r') as myfile:
        if TYPE == ("SU2_CFD_AD"or"SU2_DOT_AD"):
            data=myfile.read().replace('MATH_PROBLEM= DIRECT', 'MATH_PROBLEM= DISCRETE_ADJOINT')
        elif TYPE == ("SU2_DEF"):
            pass
        else:
            pass


def ReadUserInput(name):
    IN = {}
    infile = open(name, 'r')
    for line in infile:
      words = re.split('=| |\n|,|[|]',line)
      if not any(words[0] in s for s in ['\n', '%', ' ', '#']):
        words = list(filter(None, words))
        for i in range(0, len(words)):
            try:
                words[i] = float(words[i])
            except:
                words[i] = words[i]
        if len(words[1::1]) == 1 and isinstance(words[1], float):
            IN[words[0]] = [words[1]]
        elif len(words[1::1]) == 1 and isinstance(words[1], str):
            IN[words[0]] = words[1]
        else:
            IN[words[0]] = words[1::1]
    IN['Config_Path'] = name
    return IN

def WriteBladeConfigFile(name,IN):
    for key in IN:
        input = str(IN[key])
        output = input.replace('[','')
        output1 = output.replace(']', '')
        name.write("%s=%s\n"%(key,output1))
