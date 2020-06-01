#Change log 26/10/2019


## New features and behavior

- Added extensive documentation to all classes and functions (Numpy-style documentation)

- The units of angular design variables were changed from radians to degrees in the configuration file (the units are converted radians within the code to perform the computations)

- Added a new 2D parametrization based on camberline and thickness
	- The type of parametrization is selected in the .cfg file (see configuration file changes)
	- The names of the design variables and function calls are handled automatically within the code
	- The regression tests were generalized to use any type of 2D parametrization

- Added new functions to compute the coordinates of the hub and shroud surfaces
	These surfaces are extended to the inlet and outlet region following the start/end slopes (G1 continuity)
	 THis generalizes the TurboGrid mesh generation files for axial, radial or mixed flow


- Added interactive plots
	The interactive plot is updated when the configuration file is changed
	This is useful to visualize the geometry changes when a design variable is updated
	This is useful to set a good initial guess for the blade matching
	
- Blade matching new features
	- Interactive plot to give an initial guess
	- The matching figures are updated each iteration (display the optimization progress visually)
	- Automation of the blade matching calling the (u,v) matching each N-iterations of the DV-matching (minimize user intervention)
	


## Directory structure and new files
- Now there are 4 executables in the `bin` directory
	1. MakeBlade.py
	Reads a configuration file, generates the corresponding blade geometry and prints files with:
	a) blade surface coordinates, b) sensitivity with respect to design variables, and c) mesh generation input
	
	2. MatchBlade.py (I decided to keep the terminology blade matching after all)
	Reads a configuration file and finds the set of design variables / (u,v) parametrization that match a prescribed blade
	
	3. PlotBlade.py
	Reads a configuration file and plots the geometry of a blade in Matplotlib or Tecplot
	A new option was added to plot 3D geometries as a suface plot in Matplotlib
	A new option was added to have 2D and 3D interactive plots that are updated when the input .cfg file is modified
	
	4. OptimizeBlade.py
	This executable is left unchanged and we will need to update it to make it compatible with the new changes in ParaBlade
	
- The  `common` directory is unchanged or contains very minor changes

- The `docs` directory is unchanged or contains very minor changes. I expect to add contents to this directory soon

- The directory `TestCases` was renamed to `projects` (this is consistend with the single-word, all lowercase names of the other directories)

-  There are mayor changes in the `src` directory but they do not affect the behavior of the executables in the `bin` directory.

	**Main files** (used by the executables in the `bin` directory

	- blade_3D.py - Contains the Blade3D class that generates the blade geometry
	- blade_match.py - Contains the BladeMatch class the matches an existing geometry
	- blade_output.py - Contains the class BladeOutput that prints the coordinates, sensitivities and mesh generation files
	- blade_plot.py - Contains the class BladePlot used to visualize the generated blades


	**Secondary files** (used by Blade3D to generate the blade geometry):
	
	- blade_2D_connecting_arcs.py - Contains the Blade2DConnectingArcs class with the **old parametrization** based on connecting arcs that was already implemented in ParaBlade
	- blade_2D_camber_thickness.py - Contains the Blade2DCamberThickness class with the **new parametrization** based on camberline and thickness
	- arc_length_integration.py
	- bicubic_interpolation.py
	- bilinear_interpolation.py
	- transfinite_interpolation.py
	- bspline_curve.py
	- nurbs_curve.py
	- bspline_surface.py (not used)
	- nurbs_surface.py (not used)
	
	In addition 3 new subdirectories were added to `src`:
	
	1. `src_developer_demos` contains python scripts to test each of the secondary src files (these can be used debug or as demonstration scripts)
	2. `src_matlab_plots` contains Matlab functions to generate publication-quality 3D figures of the blade geometry (the Matplotlib 3D renderer is not as good as the Matlab one)
	3. `src_deprecated` contains pieces of deprecated code that might be recycled in the future (this directory can be ignored)
	



## Changes in the configuration file

- The configuration file was prettified with sections

- New comments were added to explain the behavior of the different options

- New notation for some of the options
	- NDIM instead of N_DIM
	- N_BLADES instead of NBLADES
	- CASCADE_TYPE instead of REP_TYPE.
	The CASCADE_TYPE option value AXISYM was renamed to ANNULAR (I think that this is more intuitive) 

- New options
	- New option N_SECTIONS. Number of sections used to interpolate the blade coordinates (this was hardcoded in Blade3D in the previous version)
	- New option PARAMETRIZATION_TYPE. This is used to specify whether to use the connecting/arcs or the camber/thickness 2D parametrization
	- New design variable names corresponding to the camber/thickness 2D parametrization
	- The option OPERATION_TYPE was reworked. New valid options: GEOMETRY  or SENSITIVITY
		GEOMETRY instructs Blade3D to compute only the surface coordinates
		SENSITIVITY instrcuts Blade3D to compute the surface coordinates, normals, and sensitivities
		(other OPERATION_TYPE options are not necessary because blade matching is done executing MatchBlade.py and blade plotting is done executing PlotBlade.py)
	


## Additiona changes
I tried to summarize the most important changes in this document but there are many other minor changes to the code. I did not document them because I think that they do not affect the "visible" behavior of ParaBlade.

Let me know if there is something that is not behaving as expected that was not documented in this code. We can also organize a "code tour" sharing screen on Skype if you want


Roberto
Saturday, 26. October 2019 05:28PM 




# Changelog 09/12/2019
- Merge of camber/thickness parametrization into master
- Changed variable name N_DIMENSIONS to NDIM
- Added interactive plot (bugs fixed)
- Updated match_blade.py after merge
- Added all the changes since 28/10/2019 that were not included during merge
- Updated regression tests after merge
- Updated plotting functions and output functions after merge
- Updated MATLAB plots after merge

Monday, 09. December 2019 07:49PM 

