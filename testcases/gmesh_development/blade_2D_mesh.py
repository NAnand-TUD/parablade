# -------------------------------------------------------------------------------------------------------------------- #
# Import general packages
# -------------------------------------------------------------------------------------------------------------------- #
import os
import sys
import pdb
import time
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gmsh_api
import gmsh_api.gmsh as gmsh


#----------------------------------------------------------------------------------------------------------------------#
# Import user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
BLADE_HOME = os.environ["BLADE_HOME"]
sys.path.append(BLADE_HOME+'/src/')
sys.path.append(BLADE_HOME+'/common/')


# -------------------------------------------------------------------------------------------------------------------- #
# Define the 2D blade class
# -------------------------------------------------------------------------------------------------------------------- #
class BladeMesh:

    """ Create a 2D blade section object

    The parametrization is based on two arcs (upper surface and lower surface)
    The arcs are constructed imposing a thickness distribution normal to the camber line
    The thickness distributions of the upper and lower surfaces are independent
    The connection between upper and lower surfaces at the leading/trailing edges is G2 continuous

    Parameters
    ----------
    blade_geo : class
    Instance of BladeGeo containing the geometry of the 2D blade and the boundaries of the flow domain

    """

    def __init__(self, blade_geo):

        # Set the input geometry as instance variable
        self.blade_geo = blade_geo


    # -------------------------------------------------------------------------------------------------------------------- #
    # Create the mesh using Gmsh
    # -------------------------------------------------------------------------------------------------------------------- #
    def make_mesh(self):


        # Initialize Gmsh. Optional inputs `argv` and `readConfigFiles`
        gmsh.initialize()

        # Clear all loaded models and post-processing data, and add a new empty model.
        gmsh.clear()

        # # If sys.argv is passed, Gmsh will parse the commandline in the same way as the standalone Gmsh app.
        # gmsh.initialize(sys.argv)

        # By default Gmsh will not print out any messages: in order to output messages
        # on the terminal, just set the standard Gmsh option "General.Terminal"
        # (same format and meaning as in .geo files):
        gmsh.option.setNumber("General.Terminal", 1)

        # Add a new model and set it as the current model
        gmsh.model.add("blade_mesh")


        # The Python API provides direct access to the internal CAD kernels
        # To create geometrical points with the built-in CAD kernel use gmsh.model.geo.addPoint():
        # - the first 3 arguments are the point coordinates (x, y, z)
        # - the next (optional) argument is the target mesh size close to the point
        # - the last (optional) argument is the point tag

        # Set a targe size for the mesh elements
        # There are many other more sophisticated ways to set mesh size using gmsh!
        lc = 3e-2

        # ---------------------------------------------- #
        # Define the flow domain
        # ---------------------------------------------- #
        #  Create the upper periodic boundary line
        upper_periodic_points = []
        Px, Py = np.real(self.blade_geo.upper_periodic_boundary.P)
        for x, y in zip(Px, Py):
            upper_periodic_points.append(gmsh.model.geo.addPoint(x, y, 0.00, lc))
        upper_periodic_line = gmsh.model.geo.addBSpline(upper_periodic_points)

        # Create the lower periodic boundary line
        lower_periodic_points = []
        Px, Py = np.real(self.blade_geo.lower_periodic_boundary.P)
        for x, y in zip(Px, Py):
            lower_periodic_points.append(gmsh.model.geo.addPoint(x, y, 0.00, lc))
        lower_periodic_line = gmsh.model.geo.addBSpline(lower_periodic_points)

        # Create the inflow/outflow boundary lines
        # The points that define inflow/outflow boundaries must be the same topological points as the periodic boundaries
        inflow_line = gmsh.model.geo.addLine(lower_periodic_points[0], upper_periodic_points[0])
        outflow_line = gmsh.model.geo.addLine(lower_periodic_points[-1], upper_periodic_points[-1])


        # ---------------------------------------------- #
        # Define the blade
        # ---------------------------------------------- #
        # Create the blade upper side line
        upper_side_points = []
        Px, Py = np.real(self.blade_geo.upper_side.P)
        for x, y in zip(Px, Py):
            upper_side_points.append(gmsh.model.geo.addPoint(x, y, 0.00, lc))
        upper_side_line = gmsh.model.geo.addBSpline(upper_side_points)

        # Create the blade lower side line
        lower_side_points = []
        Px, Py = np.real(self.blade_geo.lower_side.P)
        for x, y in zip(Px, Py):
            lower_side_points.append(gmsh.model.geo.addPoint(x, y, 0.00, lc))

        # The first and last points of the upper and lower sides of the blade must be the same topological point
        lower_side_points[0] = upper_side_points[0]
        lower_side_points[-1] = upper_side_points[-1]
        lower_side_line = gmsh.model.geo.addBSpline(lower_side_points)


        # Combine the lines into curve loops that are used to create the surface
        curve_loop_1 = gmsh.model.geo.addCurveLoop([lower_periodic_line, outflow_line, -upper_periodic_line, -inflow_line])
        curve_loop_2 = gmsh.model.geo.addCurveLoop([lower_side_line, -upper_side_line])

        # The domain is defined by the periodic/inflow/outflow curve loop with a hole defined by the blade curve loop
        domain = gmsh.model.geo.addPlaneSurface([curve_loop_1, curve_loop_2])
        # domain = gmsh.model.geo.addPlaneSurface([curve_loop_1])

        # ---------------------------------------------- #
        # Set physical groups
        # ---------------------------------------------- #
        # Physical groups are defined by providing the dimension of the group (0 for
        # physical points, 1 for physical curves, 2 for physical surfaces and 3 for
        # phsyical volumes) followed by a vector of entity tags. The last (optional)
        # argument is the tag of the new group to create.

        # Set a tag for the blade boundary curves
        blade_boundary_tag = 1
        gmsh.model.addPhysicalGroup(1, [upper_side_line, upper_side_line], blade_boundary_tag)
        gmsh.model.setPhysicalName(1, blade_boundary_tag, "Blade boundary")

        # Set a tag for the flow domain surface
        gmsh.model.addPhysicalGroup(2, [domain], 1)
        gmsh.model.setPhysicalName(2, 1, "Flow domain")


        # ---------------------------------------------- #
        # Syncronize CAD kernel
        # ---------------------------------------------- #
        # Before it can be meshed, the internal CAD representation must be synchronized
        # with the Gmsh model, which will create the relevant Gmsh data structures. This
        # is achieved by the gmsh.model.geo.synchronize() API call for the built-in CAD
        # kernel. Synchronizations can be called at any time, but they involve a non
        # trivial amount of processing; so while you could synchronize the internal CAD
        # data after every CAD command, it is usually better to minimize the number of
        # synchronization points.
        gmsh.model.geo.synchronize()


        # ---------------------------------------------- #
        # Set meshing algorithm and options
        # ---------------------------------------------- #
        # Apply an elliptic smoother to the grid
        gmsh.option.setNumber("Mesh.Algorithm", 1)
        # gmsh.option.setNumber("Mesh.Smoothing", 100)

        # Generate a mesh of the current model
        # dim possible values (0, 1, 2, or 3)
        gmsh.model.mesh.generate(dim=2)

        # Refine the mesh
        # gmsh.model.mesh.refine()

        # Recombine the triangles into quads
        # gmsh.model.mesh.recombine()
        # gmsh.model.mesh.splitQuadrangles(quality=0.7)

        # # Set the mesh optimization algorithm (3D only?)
        # # method possible values ('', 'Netgen', 'HighOrder', HighOrderElastic', 'Laplace2D')
        # gmsh.model.mesh.optimize(method='', force=False, niter=1)


        # ---------------------------------------------- #
        # Save the mesh
        # ---------------------------------------------- #
        gmsh.write("./mesh_files/blade_mesh.su2")
        gmsh.write("./mesh_files/blade_mesh.msh")
        gmsh.write("./mesh_files/blade_mesh.geo_unrolled")

        # Remember that by default, if physical groups are defined, Gmsh will export in
        # the output mesh file only those elements that belong to at least one physical
        # group. To force Gmsh to save all elements, you can use
        # gmsh.option.setNumber("Mesh.SaveAll", 1)

        # This should be called at the end:
        gmsh.finalize()


        # ---------------------------------------------- #
        # Futher development
        # ---------------------------------------------- #
        # Create a multiblock structured mesh using transfinite surfaces
        # Create inflation layers close to the wall
        #     1) Check boundary layer feature
        #     2) Add a mesh element size field with the distance to the wall
        # Explore the mesh generation algorithms and options
        # Do the blade meshing using the Open Cascade CAD kernel rather than the default Gmsh CAD kernel
        # Try the outout mesh on SU2
        # Try to deform the output mesh, not only the blade surface but also displace the periodic surfaces when the
        # blade spacing is changes (change in the number of blades)
        # Extend the methodology to 3D, big challenge!


        # gmsh.model.remove()                                # Remove current model
        # gmsh.model.list()                                  # List the names of all models
        # gmsh.model.setCurrent(name)                        # Set the current model to the model with name name
        # gmsh.model.getEntities()                           # Get all the entities in the current model
        # gmsh.model.setEntityName(dim, tag, name)           # Set the name of the entity of dimension tag
        # gmsh.model.getEntityName(dim, tag)                 # Get the name of the entity of dimension and tag
        # gmsh.model.getPhysicalGroups(dim)                  # Get all the physical groups in the current model.
        # gmsh.model.getEntitiesForPhysicalGroup(dim, tag)   # Get the tags of the model entities making up the physical
        # gmsh.model.getPhysicalGroupsForEntity(dim, tag)    # Get the tags of the physical groups (if any) to which the model entity belongs.
        # gmsh.model.addPhysicalGroup(dim,tags,tag)          # Add a physical group of grouping the model entities and set a tag
        # gmsh.model.getPhysicalName(dim, tag)               # Get the name of the physical group of dimension dim and tag tag.
        # gmsh.model.setPhysicalName(dim, tag, name)         # Set the name of the physical group of dimension dim and tag tag.








