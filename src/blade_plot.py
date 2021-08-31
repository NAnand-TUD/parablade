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

#----------------------------------------------------------------------------------------------------------------------#
# Importing general packages
#----------------------------------------------------------------------------------------------------------------------#
import os
import sys
import pdb
import time
import copy
import numpy as np


#---------------------------------------------------------------------------------------------#
# Setting Environment
#---------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]

#---------------------------------------------------------------------------------------------#
# Importing ParaBlade classes and functions
#---------------------------------------------------------------------------------------------#
from common.config import ReadUserInput
from src.blade_3D import Blade3D


#----------------------------------------------------------------------------------------------------------------------#
# "Cluster mode" imports
#----------------------------------------------------------------------------------------------------------------------#
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except:
    pass


#----------------------------------------------------------------------------------------------------------------------#
# Define BladePlot class
#----------------------------------------------------------------------------------------------------------------------#
class BladePlot:

    # Declare additional variables as instance variables
    fig_2D = None
    fig_3D = None
    ax_2D = None
    ax_3D = None
    points_2D = None
    points_3D = None
    sensitivity_normal = None

    def __init__(self, blade_in):

        # Declare blade_in variables as instance variables
        # Configuration options
        self.blade_in       = blade_in
        self.NDIM           = blade_in.NDIM
        self.N_SECTIONS     = blade_in.N_SECTIONS
        self.N_BLADES       = blade_in.N_BLADES
        self.PLOT_FORMAT    = blade_in.PLOT_FORMAT
        self.CASCADE_TYPE   = blade_in.CASCADE_TYPE
        self.OPERATION_TYPE = blade_in.OPERATION_TYPE
        self.IN             = blade_in.IN

        # Design variables
        self.DVs_names = blade_in.DVs_names
        self.DVs_functions = blade_in.DVs_functions
        self.DVs_control_points = blade_in.DVs_control_points

        # (u,v) parametrization
        self.u = blade_in.u
        self.v = blade_in.v
        self.Nu = blade_in.Nu                   # Number of u query points from meshgrid format
        self.Nv = blade_in.Nv                   # Number of u query points from meshgrid format
        self.N_points = blade_in.N_points

        # Blade geometry
        self.surface_coordinates = np.real(blade_in.surface_coordinates)
        self.surface_normals = np.real(blade_in.surface_normals)
        self.surface_sensitivity = blade_in.surface_sensitivity
        self.hub_coordinates = np.real(blade_in.hub_coordinates)
        self.shroud_coordinates = np.real(blade_in.shroud_coordinates)

        # Create output directory
        os.system("rm -rf output")
        if os.path.exists("output") is False:
            os.mkdir("output")


    # ---------------------------------------------------------------------------------------------------------------- #
    # Make default plots
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_plots(self):

        """ Call the Matplotlib or Tecplot 360 plotting functions """

        if self.PLOT_FORMAT == 'MATPLOTLIB':
            self.make_python_plot()
            plt.show()
        elif self.PLOT_FORMAT == 'INTERACTIVE':
            self.make_interactive_plot()
            plt.show()
        elif self.PLOT_FORMAT == 'TECPLOT':
            self.make_tecplot_plot()
        else:
            raise Exception('Choose a valid plotting option: "MATPLOTLIB" or "TECPLOT"')


    def make_python_plot(self):

        """ Call the Matplotlib plotting functions """

        # Import the Matplotlib library
        try:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
        except:
            raise Exception('Matplotlib library not installed... exiting... \n')

        # Call the 2D or 3D plotting functions
        if self.NDIM == 2:
            self.make_plot_matplotlib_2D()
        elif self.NDIM == 3:
            self.make_plot_matplotlib_3D()
            # self.plot_meridional_channel()
        else:
            raise Exception('The number of dimensions must be 2 or 3')


    def make_tecplot_plot(self):

        """ Call the Tecplot 360 plotting functions """

        # Import the Tecplot library
        try:
            import tecplot as tp
        except:
            raise Exception('Tecplot 360 Python library not installed... exiting... \n')

        # Call the 2D or 3D plotting functions
        if self.NDIM == 2:
            self.make_plot_tecplot_2D()
        elif self.NDIM == 3:
            self.tecplot_3D()
        else:
            raise Exception('The number of dimensions must be 2 or 3')


    # ---------------------------------------------------------------------------------------------------------------- #
    # TECPLOT 360 plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_plot_tecplot_2D(self):

        """ Create 2D Tecplot 360 plot """

        blade = open("2D_Blade.dat",'w')
        blade.write("VARIABLES = \"X\", \"Y\"\n")
        for i in self.surface_coordinates:
            blade.write("%e, \t %e \n"%(i[0],i[1]))
        blade.close()


    def tecplot_3D(self):

        """ Create 3D Tecplot 360 plot """
        # TODO this function needs some cleaning and maybe automation in Tecplot

        if self.OPERATION_TYPE == "SENSITIVITY":
            self.sensitivity_normal = {}
            for key in self.DVs_names:
                for number in range(len(self.DVs_control_points[key])):
                    self.sensitivity_normal[key + '_' + str(number)] = np.sum(
                        self.surface_normals * self.surface_sensitivity[key + '_' + str(number)], axis=0)

        # Prepare the surface plot
        x = self.surface_coordinates[0, :]
        y = self.surface_coordinates[1, :]
        z = self.surface_coordinates[2, :]

        # Prepare the surface plot
        import tecplot as tp
        from tecplot.session import set_style
        from tecplot.constant import ColorMapDistribution
        import logging as log
        header = ['x', 'y', 'z']
        if self.OPERATION_TYPE == "SENSITIVITY":
            for key in self.DVs_names:
                for number in range(len(self.DVs_control_points[key])):
                    header = header + ['Sens_'+key+'_'+str(number)]

        with tp.session.suspend():
            log.info('creating tecplot dataset')
            ds = tp.active_frame().create_dataset('Data', [header])
            sphere_zone = ds.add_ordered_zone("Blade", (self.Nu, self.Nv))
            #pdb.set_trace()
            sphere_zone.values('x')[:] = x.ravel()
            sphere_zone.values('y')[:] = y.ravel()
            sphere_zone.values('z')[:] = z.ravel()
            if self.OPERATION_TYPE == "SENSITIVITY":
                for key in self.DVs_names:
                    for number in range(len(self.DVs_control_points[key])):
                        sphere_zone.values('Sens_'+key+'_'+str(number))[:] = self.sensitivity_normal[key + '_' + str(number)]
            log.info('setting plot type to Cart3D')
            tp.active_frame().plot_type = tp.constant.PlotType.Cartesian3D
            plot = tp.active_frame().plot()
            tp.data.save_tecplot_plt('tecplot_blade.plt')



    # ---------------------------------------------------------------------------------------------------------------- #
    # MATPLOTLIB plotting functions
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_plot_matplotlib_2D(self):

        """ Create 2D Matplotlib line plot """

        # Get blade coordinates
        x = self.surface_coordinates[0, :]
        y = self.surface_coordinates[1, :]

        # Plot the prescribed and the fitted blades
        self.fig_2D = plt.figure(figsize=(8, 6))
        self.ax_2D = self.fig_2D.add_subplot(111)
        self.ax_2D.set_xlabel('$x$ axis', fontsize=12, color='k', labelpad=12)
        self.ax_2D.set_ylabel('$y$ axis', fontsize=12, color='k', labelpad=12)
        # self.ax_2D.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        # self.ax_2D.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        for t in self.ax_2D.xaxis.get_major_ticks(): t.label.set_fontsize(12)
        for t in self.ax_2D.yaxis.get_major_ticks(): t.label.set_fontsize(12)

        # Plot prescribed coordinates
        self.points_2D, = self.ax_2D.plot(x, y)
        self.points_2D.set_marker(" ")
        self.points_2D.set_markersize(3.5)
        self.points_2D.set_markeredgewidth(0.5)
        self.points_2D.set_markeredgecolor("k")
        self.points_2D.set_markerfacecolor("w")
        self.points_2D.set_linestyle("-")
        self.points_2D.set_color("k")
        self.points_2D.set_linewidth(1.50)

        # Set the aspect ratio of the data
        self.ax_2D.set_aspect(1.0)

        # Set axes aspect ratio
        x_min, x_max = self.ax_2D.get_xlim()
        y_min, y_max = self.ax_2D.get_ylim()
        x_mid = (x_min + x_max) / 2
        y_mid = (y_min + y_max) / 2
        L = np.max((x_max - x_min, y_max - y_min)) / 2
        self.ax_2D.set_xlim([x_mid - 1.25 * L, x_mid + 1.25 * L])
        self.ax_2D.set_ylim([y_mid - 1.25 * L, y_mid + 1.25 * L])

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)


    def make_plot_matplotlib_3D(self):

        """ Create 3D Matplotlib scatter plot """

        # Prepare the plot
        self.fig_3D = plt.figure(figsize=(8, 7))
        self.ax_3D = self.fig_3D.add_subplot(111, projection='3d')
        self.ax_3D.view_init(azim=-55, elev=30)
        self.ax_3D.grid(False)
        self.ax_3D.xaxis.pane.fill = False
        self.ax_3D.yaxis.pane.fill = False
        self.ax_3D.zaxis.pane.fill = False
        self.ax_3D.xaxis.pane.set_edgecolor('k')
        self.ax_3D.yaxis.pane.set_edgecolor('k')
        self.ax_3D.zaxis.pane.set_edgecolor('k')
        self.ax_3D.xaxis.pane._alpha = 0.9
        self.ax_3D.yaxis.pane._alpha = 0.9
        self.ax_3D.zaxis.pane._alpha = 0.9
        self.ax_3D.set_xlabel('$x$ axis', fontsize=12, color='k', labelpad=16)
        self.ax_3D.set_ylabel('$y$ axis', fontsize=12, color='k', labelpad=16)
        self.ax_3D.set_zlabel('$z$ axis', fontsize=12, color='k', labelpad=16)
        # self.ax_3D.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # self.ax_3D.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # self.ax_3D.zaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in self.ax_3D.xaxis.get_major_ticks(): t.label.set_fontsize(12)
        for t in self.ax_3D.yaxis.get_major_ticks(): t.label.set_fontsize(12)
        for t in self.ax_3D.zaxis.get_major_ticks(): t.label.set_fontsize(12)
        self.ax_3D.xaxis.set_rotate_label(False)
        self.ax_3D.yaxis.set_rotate_label(False)
        self.ax_3D.zaxis.set_rotate_label(False)
        # self.ax_3D.set_xticks([])
        # self.ax_3D.set_yticks([])
        # self.ax_3D.set_zticks([])
        # self.ax_3D.axis('off')

        # Blade coordinates
        x_blade = self.surface_coordinates[0, :]
        y_blade = self.surface_coordinates[1, :]
        z_blade = self.surface_coordinates[2, :]

        # # Plot blade coordinates (scatter plot)
        # if scatter == 'yes':
        #     scatter_3D, = self.ax_3D.plot(x_blade, y_blade, z_blade)
        #     scatter_3D.set_marker("o")
        #     scatter_3D.set_markersize(2.5)
        #     scatter_3D.set_markeredgewidth(0.5)
        #     scatter_3D.set_markeredgecolor("k")
        #     scatter_3D.set_markerfacecolor("w")
        #     scatter_3D.set_linestyle(" ")
        #     scatter_3D.set_color("k")
        #     scatter_3D.set_linewidth(0.50)

        # Plot the hub and shroud sections
        self.blade_sections = [None, None]
        for k, v in enumerate([0, 1]):
            u_plot = np.linspace(0, 1, 250)
            v_plot = v*np.ones((250,))
            x, y, z = np.real(self.blade_in.get_surface_coordinates(u_plot, v_plot))
            self.blade_sections[k], = self.ax_3D.plot(x, y, z)
            self.blade_sections[k].set_marker(" ")
            self.blade_sections[k].set_markersize(3.5)
            self.blade_sections[k].set_markeredgewidth(0.5)
            self.blade_sections[k].set_markeredgecolor("k")
            self.blade_sections[k].set_markerfacecolor("w")
            self.blade_sections[k].set_linestyle("-")
            self.blade_sections[k].set_color("k")
            self.blade_sections[k].set_linewidth(1.50)

        # Reshape to meshgrid format
        x = x_blade.reshape((self.Nv, self.Nu))
        y = y_blade.reshape((self.Nv, self.Nu))
        z = z_blade.reshape((self.Nv, self.Nu))

        # Plot the blade surface
        self.blade_surface = self.ax_3D.plot_surface(x, y, z,
                                                     color='blue',
                                                     # edgecolor='k',
                                                     linewidth=0.50,
                                                     alpha=0.6,
                                                     shade=False,
                                                     antialiased=True,
                                                     zorder=3)

        # Set axes aspect ratio
        x_min, x_max = self.ax_3D.get_xlim()
        y_min, y_max = self.ax_3D.get_ylim()
        z_min, z_max = self.ax_3D.get_zlim()
        x_mid = (x_min + x_max) / 2
        y_mid = (y_min + y_max) / 2
        z_mid = (z_min + z_max) / 2
        L = np.max((x_max - x_min, y_max - y_min, z_max - z_min)) / 2
        self.ax_3D.set_xlim3d(x_mid - 1.0 * L, x_mid + 1.0 * L)
        self.ax_3D.set_ylim3d(y_mid - 1.0 * L, y_mid + 1.0 * L)
        self.ax_3D.set_zlim3d(z_mid - 1.0 * L, z_mid + 1.0 * L)

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)


    # ---------------------------------------------------------------------------------------------------------------- #
    # Plot the blade in interactive mode
    # ---------------------------------------------------------------------------------------------------------------- #
    def make_interactive_plot(self):

        """ Create an iteractive plot that updates the geometry in real time when the .cfg file is changed """
        try:
            import matplotlib.animation as animation
            import matplotlib.pyplot as plt
        except:
            raise Exception('Matplotlib library not installed... exiting... \n')

        print("\n")
        print("Opening the interactive plot...\n")
        print("\tEdit the .cfg file to modify the geometry in real time")
        print("\tClose the figures when you are done to continue execution\n")

        # This calls the function animate with a refresh time of "interval" miliseconds
        if self.NDIM == 2:
            self.make_plot_matplotlib_2D()
            animation_2D = animation.FuncAnimation(self.fig_2D, self.animate, interval=500)
        elif self.NDIM == 3:
            self.make_plot_matplotlib_3D()
            animation_3D = animation.FuncAnimation(self.fig_3D, self.animate, interval=1500)
        else:
            raise Exception('The number of dimensions must be "2" or "3"')

        plt.show()
        print("Closing the interactive plot...\n")


    def animate(self, _):

        """ Update the geometry of the interactive plot """

        # Read the .cfg file and update the geometry
        while True:
            try:
                self.IN = ReadUserInput(self.IN['Config_Path'])
                self.blade_in = Blade3D(self.IN)
                self.blade_in.make_blade()
                self.surface_coordinates = self.blade_in.get_surface_coordinates(self.u, self.v)
                break
            except:
                print("The .cfg file could not be loaded, trying again...")
                time.sleep(1.5)     # Wait a moment in case the .cfg file is not found

        # Rename coordinates
        x_blade = np.real(self.surface_coordinates[0, :])
        y_blade = np.real(self.surface_coordinates[1, :])
        z_blade = np.real(self.surface_coordinates[2, :])

        # Reshape to meshgrid format
        x = x_blade.reshape((self.Nv, self.Nu))
        y = y_blade.reshape((self.Nv, self.Nu))
        z = z_blade.reshape((self.Nv, self.Nu))

        # Update the coordinate values
        if self.NDIM == 2:
            self.points_2D.set_xdata(x_blade)
            self.points_2D.set_ydata(y_blade)
        elif self.NDIM == 3:

            # Remove old plot (Matplotlib does not allow to update 3D plots)
            self.blade_surface.remove()
            for section in self.blade_sections:
                section.remove()

            # Plot the hub and shroud sections
            self.blade_sections = [None, None]
            for k, v in enumerate([0, 1]):
                u_plot = np.linspace(0, 1, 250)
                v_plot = v * np.ones((250,))
                x, y, z = np.real(self.blade_in.get_surface_coordinates(u_plot, v_plot))
                self.blade_sections[k], = self.ax_3D.plot(x, y, z)
                self.blade_sections[k].set_marker(" ")
                self.blade_sections[k].set_markersize(3.5)
                self.blade_sections[k].set_markeredgewidth(0.5)
                self.blade_sections[k].set_markeredgecolor("k")
                self.blade_sections[k].set_markerfacecolor("w")
                self.blade_sections[k].set_linestyle("-")
                self.blade_sections[k].set_color("k")
                self.blade_sections[k].set_linewidth(1.50)

            # Reshape to meshgrid format
            x = x_blade.reshape((self.Nv, self.Nu))
            y = y_blade.reshape((self.Nv, self.Nu))
            z = z_blade.reshape((self.Nv, self.Nu))

            # Plot the blade surface
            self.blade_surface = self.ax_3D.plot_surface(x, y, z,
                                                         color='blue',
                                                         # edgecolor='k',
                                                         linewidth=0.50,
                                                         alpha=0.6,
                                                         shade=False,
                                                         antialiased=True,
                                                         zorder=3)


            # self.blade_surface.set_3d_properties(z_blade)
        else:
            raise Exception('The number of dimensions must be "2" or "3"')



    # ---------------------------------------------------------------------------------------------------------------- #
    # Additional plots
    # ---------------------------------------------------------------------------------------------------------------- #
    def plot_meridional_channel(self):

        # Prepare the figure
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        fontsize = 10
        ax.set_xlabel('$z$ - axis', fontsize=fontsize, color='k', labelpad=12)
        ax.set_ylabel('$x$ - axis', fontsize=fontsize, color='k', labelpad=12)
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)

        # Plot extended hub
        x, z = self.hub_coordinates
        line, = ax.plot(x, z)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        # Plot extended shroud
        x, z = self.shroud_coordinates
        line, = ax.plot(x, z)
        line.set_linewidth(1.25)
        line.set_linestyle("-")
        line.set_color("k")
        line.set_marker(" ")
        line.set_markersize(3.5)
        line.set_markeredgewidth(1)
        line.set_markeredgecolor("k")
        line.set_markerfacecolor("w")
        line.set_label(' ')

        # Plot the blade
        u = np.linspace(0.00, 1.00, 1000)
        names = [('x_hub', 'z_hub'), ('x_shroud', 'z_shroud'), ('x_leading', 'z_leading'), ('x_trailing', 'z_trailing')]
        for name in names:

            # Plot the meridional channel boundaries
            x = np.real(self.DVs_functions[name[0]](u)[0,:])
            z = np.real(self.DVs_functions[name[1]](u)[0,:])
            line, = ax.plot(x, z)
            line.set_linewidth(1.25)
            line.set_linestyle("-")
            line.set_color("k")
            line.set_marker(" ")
            line.set_markersize(3.5)
            line.set_markeredgewidth(1)
            line.set_markeredgecolor("k")
            line.set_markerfacecolor("w")
            line.set_label(' ')

            # Plot the meridional channel control points
            Px = np.asarray(self.DVs_control_points[name[0]])
            Pz = np.asarray(self.DVs_control_points[name[1]])
            points, = ax.plot(Px, Pz)
            points.set_linewidth(0.75)
            points.set_linestyle("-.")
            points.set_color("k")
            points.set_marker("o")
            points.set_markersize(4)
            points.set_markeredgewidth(1)
            points.set_markeredgecolor("r")
            points.set_markerfacecolor("w")
            line.set_label(' ')

        # Set the aspect ratio of the data
        ax.set_aspect(1.0)

        # Adjust pad
        plt.tight_layout(pad=5.0, w_pad=None, h_pad=None)

        # Hide axes
        plt.axis('off')
