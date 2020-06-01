This directory contains regression tests to assess:

	1) The accuracy of the normal vector computation
		- Forward finite differences
		- Central finite differences
		- Complex step

	2) The accuracy of the sensitivity computation
		- Forward finite differences
		- Central finite differences
		- Complex step

	3) The continuity of the blade parametrization at the edges
		- Surface coordinates continuity (G0)
		- Normal vector continuity (G1)
		- Surface sensitivity continuity


In addition, the directory "influence_design_variables_2D" contains an script to visualize the influence of the design variables that define the 2D parametrization of the blade sections


29/05/2019
- Roberto Agromayor


The regression test directory scripts were updated to deal both with

    1) CAMBERLINE_THICKNESS parametrizations
    2) CONNECTING_ARC parametrizations

13/10/2019
- Roberto Agromayor
