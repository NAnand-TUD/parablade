# def make_spline_surface(self):
#     # This code may be kept as the bspline approximation option but has to be cleaned
#     # Get the NURBS variation of the design variables (create a NURBS function for each design variable)
#     self.get_BSpline_variation()
#
#     # Get the characteristic length (mean-line arc-length) to scale the y-coordinate
#     self.characteristic_length = self.get_characteristic_length()
#
#     h = 1e-6
#     Nu, Nv = 100, 11
#     u = np.linspace(0.00 + h, 1.00 - h, Nu)
#     v = np.linspace(0.00 + h, 1.00 - h, Nv)
#     [uu, vv] = np.meshgrid(u, v, indexing='xy')
#     self.u = uu.flatten()
#     self.v = vv.flatten()
#
#     # Array of control points
#     P = np.zeros((3, Nu, Nv), dtype=complex)
#     for k in range(Nv):
#         P[:, :, k] = self.get_section_coordinates(u, v[k])
#
#     # Maximum index of the control points (counting from zero)
#     n = np.shape(P)[1] - 1
#     m = np.shape(P)[2] - 1
#
#     # Define the order of the basis polynomials
#     # Linear (p = 1), Quadratic (p = 2), Cubic (p = 3), etc.
#     # Set p = n (number of control points minus one) to obtain a Bezier
#     p = 3
#     q = 3
#
#     # Definition of the knot vectors (clamped spline)
#     # p+1 zeros, n minus p equispaced points between 0 and 1, and p+1 ones.  In total r+1 points where r=n+p+1
#     # q+1 zeros, m minus q equispaced points between 0 and 1, and q+1 ones. In total s+1 points where s=m+q+1
#     U = np.concatenate((np.zeros(p), np.linspace(0, 1, n - p + 2), np.ones(p)))
#     V = np.concatenate((np.zeros(q), np.linspace(0, 1, m - q + 2), np.ones(q)))
#
#     # Create, evaluate, and plot the spline surface
#     my_BSpline_surface = BSplineSurface(P, p, q, U, V)
#
#     h = 1e-6
#     Nu, Nv = 1000, 500
#     u = np.linspace(0.00 + h, 1.00 - h, Nu)
#     v = np.linspace(0.00 + h, 1.00 - h, Nv)
#     [uu, vv] = np.meshgrid(u, v, indexing='xy')
#     self.u = uu.flatten()
#     self.v = vv.flatten()
#
#     surface_coordinates = my_BSpline_surface.get_BSpline_surface_value(self.u, self.v)
#
#     return surface_coordinates




# # Plot the hub and shroud surface lines in the (x,z) plane
# fig, ax = plt.subplots()
# fontsize = 12
# ax.set_xmargin(0.1)
# ax.set_ymargin(0.1)
# ax.set_aspect(1.00)
# ax.set_xlabel('z -- axis', fontsize=fontsize, color='k', labelpad=12)
# ax.set_ylabel('x -- axis', fontsize=fontsize, color='k', labelpad=12)
# for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
# for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(fontsize)
# points1, = ax.plot(self.hub_coordinates[:, 1], self.hub_coordinates[:, 0])
# points1.set_marker(" ")
# points1.set_markersize(3.5)
# points1.set_markeredgewidth(0.5)
# points1.set_markeredgecolor("k")
# points1.set_markerfacecolor("w")
# points1.set_linestyle("-")
# points1.set_color("k")
# points1.set_linewidth(0.75)
#
# points1, = ax.plot(self.shroud_coordinates[:, 1], self.shroud_coordinates[:, 0])
# points1.set_marker(" ")
# points1.set_markersize(3.5)
# points1.set_markeredgewidth(0.5)
# points1.set_markeredgecolor("k")
# points1.set_markerfacecolor("w")
# points1.set_linestyle("-")
# points1.set_color("k")
# points1.set_linewidth(0.75)
#
# plt.show()




# def plot_surface_coordinates(self, ax):
#     C = np.real(self.surface_coordinates)
#     points, = ax.plot(C[:, 0], C[:, 1], C[:, 2], label='parametric curve')
#     points.set_linestyle(" ")
#     points.set_color("k")
#     points.set_marker("o")
#     points.set_markersize(1)
#     points.set_markeredgewidth(1)
#     points.set_markeredgecolor("k")
#     points.set_markerfacecolor("w")
#
# def plot_surface_normals(self, ax):
#     size = 0.05
#     C = np.real(self.surface_coordinates)
#     N = np.real(self.surface_normals)
#     for i in range(self.N_points):
#         normals, = ax.plot(C[i, 0] + [0, N[i, 0]*size], C[i, 1] + [0, N[i, 1]*size], C[i, 2] + [0, N[i, 2]*size], label='parametric curve')
#         normals.set_linestyle("-")
#         normals.set_color("k")
#         normals.set_marker(" ")
#
# def plot_surface_sensitivity_matplotlib(self, fig, ax, u, v):
#
#     # sensitivity_normal = {}
#     # for key in self.variable_names:
#     #     for number in range(len(self.control_points[key])):
#     #         self.sensitivity_normal[key + '_' + str(number)] = np.sum(
#     #             self.surface_normals * self.surface_sensitivity[key + '_' + str(number)], axis=1)
#
#     sensitivity_normal = {}
#     key = "stagger"
#     number = 1
#
#     sensitivity_normal[key + '_' + str(number)] = np.sum(
#         self.surface_normals * self.surface_sensitivity[key + '_' + str(number)], axis=1)
#
#     # Preallocate space
#     Nv = 5
#     Nu = round(self.N_points/Nv)
#
#     X = np.zeros((Nu, Nv))
#     Y = np.zeros((Nu, Nv))
#     Z = np.zeros((Nu, Nv))
#     G = np.zeros((Nu, Nv))
#
#     # Prepare the surface plot
#     x = np.real(self.surface_coordinates[:, 0])
#     y = np.real(self.surface_coordinates[:, 1])
#     z = np.real(self.surface_coordinates[:, 2])
#     g = np.real(sensitivity_normal[key + '_' + str(number)][:])
#
#
#     # Reshape the coordinates into an array
#     for i in range(Nv):
#         X[:, i] = x[i * Nu:(i + 1) * Nu]
#         Y[:, i] = y[i * Nu:(i + 1) * Nu]
#         Z[:, i] = z[i * Nu:(i + 1) * Nu]
#         G[:, i] = g[i * Nu:(i + 1) * Nu]
#
#     # Prepare the plot
#     ax.view_init(40, 180)
#     ax.set_aspect('equal')
#     ax.set_xmargin(0.1)
#     ax.set_ymargin(0.1)
#     ax.set_zmargin(0.1)
#     # ax.grid(False)
#     # ax.set_axis_off()
#
#     # Plot the surfaces
#     ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=True,
#                     facecolors=cm.jet(G))
#
#
# # ---------------------------------------------------------------------------------------------------------------- #
# # ---------------------------------------------------------------------------------------------------------------- #
#
#
# def plot_meridional_channel(self):
#     fig, ax = plt.subplots()
#     u = np.linspace(0, 1, 100)
#     ax.set_aspect('equal')
#     x_leading = self.NURBS_variation['x_leading'](u)
#     z_leading = self.NURBS_variation['z_leading'](u)
#     x_trailing = self.NURBS_variation['x_trailing'](u)
#     z_trailing = self.NURBS_variation['z_trailing'](u)
#     x_hub = self.NURBS_variation['x_hub'](u)
#     z_hub = self.NURBS_variation['z_hub'](u)
#     x_shroud = self.NURBS_variation['x_shroud'](u)
#     z_shroud = self.NURBS_variation['z_shroud'](u)
#     ax.plot(x_leading, z_leading, 'k')
#     ax.plot(x_trailing, z_trailing, 'k')
#     ax.plot(x_hub, z_hub, 'k')
#     ax.plot(x_shroud, z_shroud, 'k')
#
#
# def plot_blade_sections(self):
#     section_variables = {}
#     fig, ax = plt.subplots()
#     for v in np.linspace(0, 1, 5):
#         # Get the design variables of the current blade section
#         for k in self.variable_names:
#             section_variables[k] = self.NURBS_variation[k].get_NURBS_value(v)
#         section = Blade2D(section_variables)
#         section.make_blade_section()
#         section.plot_blade_section(ax)
#
# def plot_NURBS_variation(self):
#     section_variables = {}
#     my_axes = {}
#     for k in self.variable_names:
#         fig, my_axes[k] = plt.subplots()
#         section_variables[k] = self.NURBS_variation[k].get_NURBS_value(np.linspace(0, 1, 100))
#         self.NURBS_variation[k].plot_NURBS_1D(ax=my_axes[k], name=k)


# if self.NDIM == 2:
#     rand_u = -0.10 + 0.20 * np.random.random(5)
#     # my_u0 = [0.100+rand_u[0], 0.300+rand_u[1], 0.500+rand_u[2], 0.700+rand_u[3], 0.900+rand_u[4]]
#     my_u0 = [0.100, 0.250, 0.500, 0.750, 0.900]
#     my_v0 = [0.500, 0.500, 0.500, 0.500, 0.500]
#     my_bounds = [[(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)]]
# else:
#     rand_u = -0.10 + 0.20*np.random.random(10)
#     rand_v = -0.25 + 0.50*np.random.random(10)
#     my_u0 = [0.100+rand_u[0], 0.300+rand_u[1], 0.500+rand_u[2], 0.700+rand_u[3], 0.900+rand_u[4],
#              0.100+rand_u[5], 0.300+rand_u[6], 0.500+rand_u[7], 0.700+rand_u[8], 0.900+rand_u[9]]
#     my_v0 = [0.250+rand_v[0], 0.250+rand_v[1], 0.250+rand_v[2], 0.250+rand_v[3], 0.250+rand_v[4],
#              0.750+rand_v[5], 0.750+rand_v[6], 0.750+rand_v[7], 0.750+rand_v[8], 0.750+rand_v[9]]
#     my_bounds = [[(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)],
#                  [(0.00 + h, 1.00 - h), (0.00 + h, 1.00 - h)]]





