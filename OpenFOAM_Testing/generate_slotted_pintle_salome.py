#%% Imports
import math
import numpy as np
import salome
from salome.geom import geomBuilder # type: ignore
from salome.smesh import smeshBuilder # type: ignore

salome.salome_init()
geompy = geomBuilder.New()
smeshpy = smeshBuilder.New()

#%% Inputs

name = "Slot_Pintle_Fluid_Domain"   # file and project name

# Fuel side
D_po = 8.0e-3   # [m], pintle post diamater
t_an = 0.5e-3   # [m], fuel annulus thickness
L_po = 10e-3    # [m], pintle post length

# Oxidizer side
D_pr = 3.0e-3   # [m], pintle rod diameter
D_mp = 4.55e-3  # [m], minimum internal diamater of pintle post
th_t = 30.0     # [°], top interal pintle angle
th_b = 20.0     # [°], bottom interal pintle angle
L_op = 1.0e-3   # [m], pintle slit opening
t_po = 0.5e-3   # [m], post thickness at slit
t_pt = 1e-3     # [m], pintle plate thickness
r_bl = 3e-3     # [m], blend radius on rod

# Domain dimensions
L_ex = 10e-3    # [m], extension of ox and f channels in -x
D_do = 70e-3   # [m], domain external wall diameter
L_do = 100e-3   # [m], domain length

# Mesh settings
lc = 1.0e-3     # [m], base mesh size
rp = 0.2        # [-], refinement factor near pintle
cf = 3.0        # [-], coarsening far away

#%% Define geometry points

tand = lambda deg: np.tan(np.deg2rad(deg))
a_bl = (np.pi/2 - np.deg2rad(th_b))/2   # blend half angle

points1 = [
    (L_po+L_op+t_pt, 0, 0),
    (L_do, 0, 0),
    (L_do, D_do/2, 0),
    (0, D_do/2, 0),
    (0, D_po/2+t_an, 0),
    (-L_ex, D_po/2+t_an, 0),
    (-L_ex, D_po/2, 0),
    (L_po, D_po/2, 0),
    (L_po, D_po/2-t_po, 0),
    (L_po - (D_po/2-t_po-D_mp/2)/tand(th_t), D_mp/2, 0),
    (-L_ex, D_mp/2, 0),
    (-L_ex, D_pr/2, 0),
    (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2, 0)]

point_bl_c = (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2+r_bl, 0)

points2 = [
    (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl)+r_bl*np.sin(2*a_bl),D_pr/2+r_bl*(1-np.cos(2*a_bl)), 0),
    (L_po+L_op, D_po/2, 0),
    (L_po+L_op+t_pt, D_po/2, 0)]

#%% Build geometry

# Define vertices
vert1 = [geompy.MakeVertex(*pt) for pt in points1]
v_center = geompy.MakeVertex(*point_bl_c)
vert2 = [geompy.MakeVertex(*pt) for pt in points2]

# Connect first part with lines
lin1 = [geompy.MakeLineTwoPnt(vert1[i], vert1[i+1]) for i in range(len(vert1) - 1)]
arc  = [geompy.MakeArcCenter(v_center, vert1[-1], vert2[0], False)]
lin2 = [geompy.MakeLineTwoPnt(vert2[i], vert2[i+1]) for i in range(len(vert2) - 1)]
closing_line = [geompy.MakeLineTwoPnt(vert2[-1], vert1[0])]

# Make closed surface
wire = geompy.MakeWire(lin1+arc+lin2+closing_line)
face = geompy.MakeFaceWires([wire], True)

#%% Revolve 90 degrees around X-axis

axis_x = geompy.MakeVectorDXDYDZ(1.0, 0.0, 0.0)
solid_90deg = geompy.MakeRevolution(face, axis_x, math.radians(90.0))
geompy.addToStudy(solid_90deg, "Fluid_Domain_3D_90deg")

#%% Find surfaces

def get_revolved_faces(solid, lines, axis, angle_rad):
    matched_faces = []
    for line in lines:
        temp_surf = geompy.MakeRevolution(line, axis, angle_rad)
        face = geompy.GetInPlace(solid, temp_surf)
        matched_faces.append(face)
    return matched_faces

# Origin lines
    # Wall lines
lin_w = lin1[2:4]
lin_w_f = [lin1[i] for i in [4,6]]
lin_w_o = [lin1[i] for i in [7,8,9,11]] + arc + lin2
    # In and outlet lines
lin_in_f = [lin1[5]]
lin_in_o = [lin1[10]]
lin_out  = [lin1[1]]

# Find corresponding faces
    # Combustion chamber faces
faces_w = get_revolved_faces(solid_90deg, lin_w, axis_x, math.radians(90))
    # Injector fuel-side faces
faces_w_f = get_revolved_faces(solid_90deg, lin_w_f, axis_x, math.radians(90))
    # Injector oxidizer-side faces
faces_w_o = get_revolved_faces(solid_90deg, lin_w_o, axis_x, math.radians(90))
    # Injector fuel inlet
faces_in_f = get_revolved_faces(solid_90deg, lin_in_f, axis_x, math.radians(90))
    # Injector oxidizer inlet
faces_in_o = get_revolved_faces(solid_90deg, lin_in_o, axis_x, math.radians(90))
    # Chamber outlet
faces_out = get_revolved_faces(solid_90deg, lin_out, axis_x, math.radians(90))

# Create geom groups
    # Combustion chamber group
grp_w = geompy.CreateGroup(solid_90deg, geompy.ShapeType["FACE"])
geompy.UnionList(grp_w, faces_w)
geompy.addToStudyInFather(solid_90deg, grp_w, "wall_chamber")
    # Injector wall group
grp_w_inj = geompy.CreateGroup(solid_90deg, geompy.ShapeType["FACE"])
geompy.UnionList(grp_w_inj, faces_w_o+faces_w_f)
geompy.addToStudyInFather(solid_90deg, grp_w_inj, "wall_injector")
    # Injector wall group
grp_out = geompy.CreateGroup(solid_90deg, geompy.ShapeType["FACE"])
geompy.UnionList(grp_out, faces_out)
geompy.addToStudyInFather(solid_90deg, grp_out, "outlet")

#%% Setup mesh

def set_mesh(type, max_size, min_size, growth_rate, group = None):
    mesh_handle = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D, geom=group) if type == "2D" else mesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    params_handle = mesh_handle.Parameters()
    params_handle.SetFineness(smeshBuilder.Custom)
    params_handle.SetMaxSize(max_size)
    params_handle.SetMinSize(min_size)
    params_handle.SetGrowthRate(growth_rate)
    params_handle.SetUseSurfaceCurvature(False)
    return mesh_handle

# Initialize mesh
mesh = smeshpy.Mesh(solid_90deg)

min_fine_size = (rp * lc) / 5.0
max_coarse_size = cf*lc

# Global element settings
mesh_3D = set_mesh("3D", max_coarse_size, min_fine_size, 0.2)
# 2D surface mesh refinement near injector
mesh_2D_inj = set_mesh("2D", rp*lc, min_fine_size, 0.15, grp_w_inj)
# 2D surface mesh coarsening on outlet
mesh_2D_out = set_mesh("2D", cf*lc, lc, 0.25, grp_out)

#%% Add inflation layers to injector and chamber walls

def add_viscous_layers(bl_thickness, no_layers, growth_ratio, face_id):
    viscous_layers = smeshpy.CreateHypothesis('ViscousLayers')
    viscous_layers.SetTotalThickness(bl_thickness)
    viscous_layers.SetNumberLayers(no_layers)
    viscous_layers.SetStretchFactor(growth_ratio)
    viscous_layers.SetFaces(face_id, 0)
    mesh.AddHypothesis(viscous_layers)
    return viscous_layers

def get_faces_id(faces: list):
    faces_id = [geompy.GetSubShapeID(solid_90deg, face) for face in faces]
    return faces_id

# Chamber wall inflation layers
visc_w = add_viscous_layers(4.03e-4, 15, 1.01, get_faces_id(faces_w))

# Injector wall inflation layers
visc_w_inj = add_viscous_layers(2.17e-5, 15, 1.13, get_faces_id(faces_w_f+faces_w_o))


#%% Compute mesh
is_done = mesh.Compute()

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()