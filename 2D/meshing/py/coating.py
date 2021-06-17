# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial for meshing a layered sediment basin stratigraphy
#
#  Geometry basics, elementary entities, physical groups
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import numpy
import gmsh

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# By default Gmsh will not print out any messages: in order to output messages
# on the terminal, just set the "General.Terminal" option to 1:
gmsh.option.setNumber("General.Terminal", 1)

# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
modelname = "coating"
gmsh.model.add(modelname)

# The Python API provides direct access to each supported geometry kernel. The
# built-in kernel is used in this first tutorial: the corresponding API
# functions have the `gmsh.model.geo' prefix.

# The first type of `elementary entity' in Gmsh is a `Point'. To create a point
# with the built-in geometry kernel, the Python API function is
# gmsh.model.geo.addPoint():
# - the first 3 arguments are the point coordinates (x, y, z)
# - the next (optional) argument is the target mesh size (the "characteristic
#   length") close to the point
# - the last (optional) argument is the point tag (a stricly positive integer
#   that uniquely identifies the point)

# The distribution of the mesh element sizes will be obtained by interpolation
# of these characteristic lengths throughout the geometry. Another method to
# specify characteristic lengths is to use general mesh size Fields (see
# `t10.py'). A particular case is the use of a background mesh (see `t7.py').
#
# If no target mesh size of provided, a default uniform coarse size will be used
# for the model, based on the overall model size.

# Units: meter
m = 0.001 	# m in mm
r0 = 1.05*m  # ball radius
r1 = 1.055*m	# coating outer radius
r2 = 100*m	# water body outer radius
#lc = 3e-1*m    # m Netzfeinheit
lc0 = 3e-1*m    # m Netzfeinheit
lc1 = 2e-3*m    # m Netzfeinheit
lc2 = 5*m    # m Netzfeinheit

# We can then define some additional points. All points should have different tags:

#gmsh.model.geo.addPoint(r1, 0, 0, lc1, 1)
#gmsh.model.geo.addPoint(r2, 0, 0, lc2, 2) # Auxiliary points for the boundary
#gmsh.model.geo.addPoint( 0,r2, 0, lc2, 3)
#gmsh.model.geo.addPoint( 0,r1, 0, lc1, 4) # Auxiliary points for the boundary

# If the tag is not provided explicitly, a new tag is automatically created, and
# returned by the function:

O = gmsh.model.geo.addPoint( 0, 0, 0, lc0)
A = gmsh.model.geo.addPoint(r0, 0, 0, lc1)
B = gmsh.model.geo.addPoint(r1, 0, 0, lc1)
C = gmsh.model.geo.addPoint(r2, 0, 0, lc2)
D = gmsh.model.geo.addPoint( 0,r2, 0, lc2)
E = gmsh.model.geo.addPoint( 0,r1, 0, lc1)
F = gmsh.model.geo.addPoint( 0,r0, 0, lc1)

# Curves are Gmsh's second type of elementery entities, and, amongst curves,
# straight lines are the simplest. The API to create straight line segments with
# the built-in kernel follows the same conventions: the first 2 arguments are
# point tags (the start and end points of the line), and the last (optional one)
# is the line tag.
#
# Note that curve tags are separate from point tags - hence we can reuse tag `1'
# for our first curve. And as a general rule, elementary entity tags in Gmsh
# have to be unique per geometrical dimension.
#gmsh.model.geo.addLine(1, 2, 1) # (start, end, tag)

gmsh.model.geo.addLine(O, A, 1)
gmsh.model.geo.addLine(A, B, 2)
gmsh.model.geo.addLine(B, C, 3)

gmsh.model.geo.addLine(D, E, 4)
gmsh.model.geo.addLine(E, F, 5)
gmsh.model.geo.addLine(F, O, 6)

gmsh.model.geo.addCircleArc(A, O, F, 7)
gmsh.model.geo.addCircleArc(B, O, E, 8)
gmsh.model.geo.addCircleArc(C, O, D, 9)

# The third elementary entity is the surface. In order to define a surface 
# from the curves defined above, a curve loop has first to be defined.
# A curve loop is defined by an ordered list of connected curves,
# a sign being associated with each curve (depending on the orientation of the
# curve to form a loop). The API function to create curve loops takes a list
# of integers as first argument, and the curve loop tag (which must be unique
# amongst curve loops) as the second (optional) argument:
# e.g. gmsh.model.geo.addCurveLoop([4, 1, 2, 3], 1)

gmsh.model.geo.addCurveLoop([1, 7, 6], 1)
gmsh.model.geo.addCurveLoop([2, 8, 5, -7], 2)
gmsh.model.geo.addCurveLoop([3, 9, 4, -8], 3)

# Add plane surfaces defined by one or more curve loops. The first curve 
# loop defines the exterior contour; additional curve loop define holes.
# (only one here, representing the external contour, since there are no holes
# --see `t4.py' for an example of a surface with a hole):
for l in range(1, 4):
	gmsh.model.geo.addPlaneSurface([l], l)

# At this level, Gmsh knows everything to display the surfaces and
# to mesh it. An optional step is needed if we want to group elementary
# geometrical entities into more meaningful groups, e.g. to define some
# mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
# material ("steel", "carbon") properties.
#
# Such groups are called "Physical Groups" in Gmsh. By default, if physical
# groups are defined, Gmsh will export in output files only mesh elements that
# belong to at least one physical group. (To force Gmsh to save all elements,
# whether they belong to physical groups or not, set the `Mesh.SaveAll' option
# to 1.) Physical groups are also identified by tags, i.e. stricly positive
# integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
# groups can also be given names.
#
# Here we define physical curves that groups the left, bottom, top and right
# curves in single groups for the one-dimensional boundary meshes

dim=0
centre = gmsh.model.addPhysicalGroup(dim, [O])
gmsh.model.setPhysicalName(dim, centre, "centre")


dim=1
bottom_ball = gmsh.model.addPhysicalGroup(dim, [1])
gmsh.model.setPhysicalName(dim, bottom_ball, "bottom_ball")

bottom_coat = gmsh.model.addPhysicalGroup(dim, [2])
gmsh.model.setPhysicalName(dim, bottom_coat, "bottom_coat")

bottom_water = gmsh.model.addPhysicalGroup(dim, [3])
gmsh.model.setPhysicalName(dim, bottom_water, "bottom_water")

left_ball = gmsh.model.addPhysicalGroup(dim, [6])
gmsh.model.setPhysicalName(dim, left_ball, "left_ball")

left_coat = gmsh.model.addPhysicalGroup(dim, [5])
gmsh.model.setPhysicalName(dim, left_coat, "left_coat")

left_water = gmsh.model.addPhysicalGroup(dim, [4])
gmsh.model.setPhysicalName(dim, left_water, "left_water")

ball = gmsh.model.addPhysicalGroup(dim, [7])
gmsh.model.setPhysicalName(dim, ball, "ball")

coat = gmsh.model.addPhysicalGroup(dim, [8])
gmsh.model.setPhysicalName(dim, coat, "coat")

water = gmsh.model.addPhysicalGroup(dim, [9])
gmsh.model.setPhysicalName(dim, water, "water")


# ... and we define physical surfaces for the different
# material domains containing the corresponding geometrical surfaces:
dim=2
Lay1 = gmsh.model.addPhysicalGroup(dim, [1])
gmsh.model.setPhysicalName(dim, Lay1, "Layer1")

Lay2 = gmsh.model.addPhysicalGroup(dim, [2])
gmsh.model.setPhysicalName(dim, Lay2, "Layer2")

Lay3 = gmsh.model.addPhysicalGroup(dim, [3])
gmsh.model.setPhysicalName(dim, Lay3, "Layer3")

# Before it can be meshed, the internal CAD representation must be synchronized
# with the Gmsh model, which will create the relevant Gmsh data structures. This
# is achieved by the gmsh.model.geo.synchronize() API call for the built-in
# geometry kernel. Synchronizations can be called at any time, but they involve
# a non trivial amount of processing; so while you could synchronize the
# internal CAD data after every CAD command, it is usually better to minimize
# the number of synchronization points.
gmsh.model.geo.synchronize()

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# ... and save it to disk: msh=current format, msh22 older format
gmsh.write("coating.msh")

# Remember that by default, if physical groups are defined, Gmsh will export in
# the output mesh file only those elements that belong to at least one physical
# group. To force Gmsh to save all elements, you can use
#
# gmsh.option.setNumber("Mesh.SaveAll", 1)

# By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
# format (the `MSH' format). You can save meshes in other mesh formats by
# specifying a filename with a different extension. For example

#gmsh.write("basinMesh.pdf") only accessible via GUI
#gmsh.write("basinMesh.pvtu")
#gmsh.write("quarterCircle.vtk")
#gmsh.write("quarterCircle.unv")
#gmsh.write("quarterCircle.stl")

# will save the mesh in the vtk, unv and stl format.

# To visualize the model we can run the graphical user interface with:
gmsh.fltk.run()


# Note that starting with Gmsh 3.0, models can be built using other geometry
# kernels than the default "built-in" kernel. To use the OpenCASCADE geometry
# kernel instead of the built-in kernel, you should use the functions with the
# `gmsh.model.occ' prefix.
#
# Different geometry kernels have different features. With OpenCASCADE, instead
# of defining the surface by successively defining 4 points, 4 curves and 1
# curve loop, one can define the rectangular surface directly with
#
# gmsh.model.occ.addRectangle(.2, 0, 0, .1, .3)
#
# Boolean operation using OpenCascade module (occ, see below)
# gmsh.model.occ.cut([(3, 1)], [(3, 2)], 3)
#
# After synchronization with the Gmsh model with
#
# gmsh.model.occ.synchronize()
#
# the underlying curves and points could be accessed with
# gmsh.model.getBoundary().
#
# See e.g. `t16.py', `t18.py', `t19.py' or `t20.py' for complete examples based
# on OpenCASCADE, and `demos/api' for more.

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
