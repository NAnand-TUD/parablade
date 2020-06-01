# ParaBlade

## SU2 Installation

If you want to know what is SU2. Please go to the link below:

[![SU2_Website](https://img.shields.io/badge/Website-SU2-blue.svg)](https://su2code.github.io/)   
[![SU2_Website](https://img.shields.io/badge/Foundation-SU2-blue.svg)](https://su2code.github.io/)

## Download SU2

To install SU2 write the following in your terminal
```
git clone https://github.com/su2code/SU2.git
```

The blade parametrization is coupled with the feature_turbomachinery branch of SU2. To change the brach use:
```
git checkout feature_turbomachinery
``` 

## Install SU2

To install SU2 make a file named ```my_preconfig_AD``` in the SU2 folder and add the following in it:
```
./preconfigure.py --prefix=<address to SU2> --enable-mpi  --with-cc=<address to mpicc>  --with-cxx=<location to mpicxx> --enable-autodiff
```

In order to make it executable, open a Terminal where the file that you just created is located and type the following command:
```
chmod +x my_precofig_AD
```

Then run ```./my_preconfig_AD``` in the Terminal and follow the instruction on the screen. Please note that the script 
`preconfigure.py` is written in Python 2.7, therefore when executing it you should be pointing to it on the header of the file.

After, run ```make -j <number of processors> install``` and wait for installation to complete.

**IMPORTANT!**

In order to be able to run SU2 using Terminal, you have to update the PATH variable. In order to do so, open your ``.bashrc`` 
(Linux) or ```.bash_profile``` (macOS) file and add the following lines:

```
export SU2_RUN="/path/to/SU2/executables"
export SU2_HOME="/path/to/SU2"
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$PYTHONPATH:$SU2_RUN
```


## Check if the installation was successful
Open a new terminal and type ```which SU2_CFD``` and if it is where you installed then the installation was successful.

## This is too difficult for me... Can I just do something more easy

Just downlaod this file and run it in your terminal.

[RunMeToGetSU2](url)

**Enjoy Optimizing Turbomachines !!!**



   
      
     
*Author: Nitish Anand, Pablo Garrido de la Serna*