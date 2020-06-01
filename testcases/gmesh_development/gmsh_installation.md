# Gmsh installation guide

1)  Download the binary files from the (official site)[http://gmsh.info/]

2) Add the path to the gmsh binaries into your `.bashrc` file.

This file is hidden in your `home` directory and you can edit it with any text editor:
		
		 cd /home/your_username/
		 gedit .bashrc

Paste the path to the gmsh binaries into the `.bashrc` file:

		export GMSH_RUN="/home/your_path_to_gmsh/gmsh-4.4.1/bin"
		export PATH=$PATH:$GMSH_RUN
		
Save `ctrl+s` and close the `.bashrc` file

To make the changes effective without closing the terminal type:
	
		source .bashrc
	


		
3) At this point you should be able to call Gmsh using the Python API (Application Programming Interface) importing the following packages at the beginning of your Python script:

		import gmsh_api
		import gmsh_api.gmsh as gmsh


4) Try to run `msh_development/main.py` to see if everything works.

