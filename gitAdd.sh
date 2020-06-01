#!/usr/bin/env sh

# Only add *.py  -- python files
#	   *.cfg -- test case configuration files
#	   *.txt -- text files for blade matching
#	   *.md  -- README file
#	   *.sh  -- shell command files

# Excluded  *.pyc -- python binary files
#	    *.csv -- output files from the code
#       *.crv -- turbogrid mesh generation files
#	    *.plt -- plotting files from tecplot
#	    *.eps -- image files
#       *.png -- image files
#       *.pdf -- image files

git add *.py *.cfg *.md *.txt *.sh *.m
