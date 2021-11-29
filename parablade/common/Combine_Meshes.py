#!/usr/bin/python3

import numpy as np
import re
import pdb

def InsertFileData(char,out):
	infile=open(char,'r')
	for line in infile:
		out.write(line)
	infile.close()

def WriteCombinedMesh(DIR,case='turb'):
	outfile=open(DIR+'/meshout_stage.su2','w')
	outfile.write('NZONE= 2\n')
	if case=='turb':
		outfile.write('IZONE= 1\n')
		InsertFileData(DIR+'/STATOR/meshout_deformed.su2',outfile)
		outfile.write('IZONE= 2\n')
		InsertFileData(DIR+'/ROTOR/meshout_deformed.su2',outfile)
	elif case=='comp':
		outfile.write('IZONE= 1\n')
		InsertFileData(DIR+'/ROTOR/meshout_deformed.su2',outfile)
		outfile.write('IZONE= 2\n')
		InsertFileData(DIR+'/STATOR/meshout_deformed.su2',outfile)
	outfile.close()
