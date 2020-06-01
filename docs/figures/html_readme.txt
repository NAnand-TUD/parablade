# ParaBlade


ParaBlade is a Python tool that implements a unified parametrization method to describe the geometry of a wide range of turbomachinery blades including:
- 2D and 3D cases
- Axial, radial and mixed flow machines
- Turbines and compressors



![python 3.2](https://img.shields.io/badge/version-latest-blue.svg)
    
![python 3.2](https://img.shields.io/badge/python-3.6.3-blue.svg)

![platform Linux,_macos,_win64](https://img.shields.io/badge/platform-Linux,_macos,_win64-blue.svg)


<div class="image123">
    <div class="imgContainer">
        <img src="./docs/figures/axial_turbine_stage.png" height="250" width="250"/>
        <p> <b> Axial flow </b> </p>
    </div>
    <div class="imgContainer">
        <img class="middle-img" src="./docs/figures/radial_outflow.png"/ height="250" width="250"/>
        <p> <b> Radial flow </b> </p>
    </div>
    <div class="imgContainer">
         <img src="./docs/figures/radial_inflow.png"/ height="" width="250"/>
        <p> <b> Mixed flow </b> </p>
    </div>
</div>




# Parametrization

The geometry of the blade is described by the parametrization of the meridional channel and the blade sections stacked along the span



## Meridional channel parametrization

The meridional channel is described by a set of four B-Splines that define 1) leading edge, 2) trailing edge, 3) hub surface and 4) shroud surface

<div class="image12">
    <div class="imgContainer">
        <img src="./docs/figures/meridional_channel_axis.svg" height="350" width="350"/>
        <p> <b> Meridional channel view </b> </p>
    </div>
    <div class="imgContainer">
        <img src="./docs/figures/meridional_channel_Bspline.svg" height="350" width="350"/>
        <p> <b> Meridional channel view </b> </p>
    </div>
</div>



## Blade section parametrization

The blade sections are defined by a set of B-Spline or NURBS segments. The control points of these curves are computed using engineering parameters such as metal angles and thickness distribution.  
There are two available parametrizations for the blade sections: 

- Connecting arcs (G1 continuous)
- Camberline and thickness (G2 continuous) 

<div class="image12">
    <div class="imgContainer">
        <img src="./docs/figures/parametrization_G1.svg" height="350" width="350"/>
        <p> <b> Connecting arcs parametrization </b> </p>
    </div>
    <div class="imgContainer">
        <img src="./docs/figures/parametrization_G2.svg" height="350" width="350"/>
        <p> <b> Camberline and thickness parametrization </b> </p>
    </div>
</div>


Both section parametrizations offer a rich design space that covers a wide range of turbomachinery blades including compressor blades as well as reaction and impulse turbine blades.

<div class="image1">
    <div class="imgContainer">
        <img src="./docs/figures/blade_morphing.gif">
    </div>
</div>


## CAD sensitivity feature

ParaBlade is able to provide the sensitivity of the surface with respect to the design variables using the complex step method. This information is required to solve shape optimization problems (e.g. maximize the blade isentropic efficiency) using gradient-based algortithms

<div class="image1">
    <div class="imgContainer">
        <img src="./docs/figures/surface_sensitivity_example.png" height="350" width="550"/>
    </div>
</div>



## Blade matching feature

ParaBlade is capable to solve the inverse problem of finding the parametrization of an existing turbomachinery blade

<div class="image12">
    <div class="imgContainer">
        <img src="./docs/figures/blade_matching_diagram.svg" height="350" width="350"/>
    </div>
    <div class="imgContainer">
        <img src="./docs/figures/blade_matching.gif" >
    </div>
</div>




# Pre-requisites

Important: MAC users, please use pip to install python packages as anaconda can give conflicts when using tecplot library.
### Pip3  
```
sudo apt-get install python-setuptools python-dev build-essential
```

### MatPlotLib
Use pip3 to install CoolProp. For more information on MatplotLib visit [here].  
[![Link matplotlib](https://img.shields.io/badge/Link-matplotlib-red.svg)](https://matplotlib.org/)
  
```
sudo pip3 install matplotlib
```
### NumPy
Use pip3 to install numpy. For more information on SciPy visit   
[![Link matplotlib](https://img.shields.io/badge/Link-numpy-red.svg)](https://www.numpy.org/).
  
```
sudo pip3 install numpy
```

### SciPy
Use pip3 to install CoolProp. For more information on SciPy visit [here]  
[![Link matplotlib](https://img.shields.io/badge/Link-scipy-red.svg)](https://www.scipy.org/).
  
```
sudo pip3 install scipy
```

### Slack Client
If you wish to get notification on slack for your optimizaiton. Please install the slack client else, just continue using
the code as it is. It would have no influence.
Use pip3 to install slack-client. For more information on slack-client visit [here].  
[![Link matplotlib](https://img.shields.io/badge/Link-slackclient-red.svg)]()
  
```
sudo pip3 install slackclient
```

Also add 
```
export SLACK_API_TOKEN="<slack legacy-token>"
```
in your ~/.bashrc file

To generate a slack legacy-tokens please visit   
[![Link matplotlib](https://img.shields.io/badge/Link-legacytoken-red.svg)](https://api.slack.com/custom-integrations/legacy-tokens)

### TecPlot
Use pip3 to install pytecplot. Information on installation is available below.  
  
```
sudo pip3 install pytecplot
```
[![Link matplotlib](https://img.shields.io/badge/Link-pytecplot-red.svg)](https://www.tecplot.com/docs/pytecplot/)

Add the tecplot libraries in your LIBRARY_PATH. Follow the instruction on pytecplot page. Link [https://www.tecplot.com/docs/pytecplot/install.html#id2]  
Important: Mac users please visit above link and read about MAC DYLD_LIBRARY_PATH
# Install
Run the command below
```
RunMe.sh
```
and follow the instructions... This should load MakeBlade.py in your working environment. 

# Usage
### TEST CASES
Go to any of the cases in TestCases Folder and run:

### [How to install SU2 ?](docs/SU2.md)

### [How to match a Blade ?](docs/Match.md)

### [How to do shape optimization ?](docs/ShapeOptimization.md)

# Branches
1. Master
2. feature_blade_matching

# Developers
The tool was developed in Power & Propulsion Group of TU Delft. The principle developers are:  
1. **Roberto** - PhD Researcher, NTNU.   
2. **Nitish Anand** - PhD Researcher, TU Delft.    

[![Link MailTo](https://img.shields.io/badge/MailTo-developers-blue.svg)](mailto:roberto.agromayor@ntnu.no;n.anand@tudelft.nl?subject=ParaBlade:Query)
# Acknowledgements
**Matteo Pini** - Assistant Prof., Power and Propulsion, TU Delft.     

# Citations
The details of the work has been published in the following documents:  

# Licence

# Standards
PEP 8
Documentation Standard: (Roberto)
  
# Some Developer guidelines

1. 

*Author: Nitish Anand*
