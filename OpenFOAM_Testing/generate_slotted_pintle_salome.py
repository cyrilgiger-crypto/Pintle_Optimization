#%% Imports
import math
import numpy as np
from pathlib import Path
from Functions.salomeToOpenFOAM import exportToFoam
import salome
from salome.geom import geomBuilder # type: ignore
from salome.smesh import smeshBuilder # type: ignore
import SMESH

salome.salome_init()
geompy = geomBuilder.New()
smeshpy = smeshBuilder.New()

#%% Inputs

name = "Slot_Pintle_Fluid_Domain"   # file and project name

# Fuel side
D_po = 8.0e-3   # [m], pintle post diamater
t_an = 0.5e-3   # [m], fuel annulus thickness
L_po = 10e-3    # [m], pintle post length
l_ch = 0.5e-3   # [m], chamfer length for fuel exit to prevent numerical instability

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
D_do = 70e-3    # [m], domain external wall diameter
L_do = 100e-3   # [m], domain length

# Mesh settings
lc = 1.0e-3     # [m], base mesh size
rp = 0.2        # [-], refinement factor near pintle
cf = 3.0        # [-], coarsening far away
angle = 10      # [°], revolution angle for periodic BC

# Misc
mesh_export = True

#%% Define geometry points

tand = lambda deg: np.tan(np.deg2rad(deg))
a_bl = (np.pi/2 - np.deg2rad(th_b))/2   # blend half angle

points1 = [
    (L_po+L_op+t_pt, 0, 0),
    (L_do, 0, 0),
    (L_do, D_do/2, 0),
    (0, D_do/2, 0),
    (0, D_po/2+t_an+l_ch, 0),
    (-l_ch, D_po/2+t_an, 0),
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

#%% Revolve x-degrees around X-axis

axis_x = geompy.MakeVectorDXDYDZ(1.0, 0.0, 0.0)
solid_rev = geompy.MakeRevolution(face, axis_x, math.radians(angle))
geompy.addToStudy(solid_rev, "Fluid_Domain_3D_rev")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()

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
lin_w_f = [lin1[i] for i in [4,5,7]]
lin_w_o = [lin1[i] for i in [8,9,10,12]] + arc + lin2
lin_w_i = closing_line
    # In and outlet lines
lin_in_f = [lin1[6]]
lin_in_o = [lin1[11]]
lin_out  = [lin1[1]]

# Find corresponding faces
    # Combustion chamber faces
faces_w = get_revolved_faces(solid_rev, lin_w, axis_x, math.radians(angle))
    # Injector fuel-side wall faces
faces_w_f = get_revolved_faces(solid_rev, lin_w_f, axis_x, math.radians(angle))
    # Injector oxidizer-side wall faces
faces_w_o = get_revolved_faces(solid_rev, lin_w_o, axis_x, math.radians(angle))
    # Injector remaining wall faces
faces_w_i = get_revolved_faces(solid_rev, lin_w_i, axis_x, math.radians(angle))
    # Injector fuel inlet
faces_in_f = get_revolved_faces(solid_rev, lin_in_f, axis_x, math.radians(angle))
    # Injector oxidizer inlet
faces_in_o = get_revolved_faces(solid_rev, lin_in_o, axis_x, math.radians(angle))
    # Chamber outlet
faces_out = get_revolved_faces(solid_rev, lin_out, axis_x, math.radians(angle))

# Find periodic BC faces at 0° and x°
faces_p_0  = geompy.GetInPlace(solid_rev, face)
faces_p_rev = geompy.GetInPlace(solid_rev, geompy.MakeRotation(face, axis_x, math.radians(angle)))

# Create geom groups
    # Combustion chamber group
grp_w = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_w, faces_w)
geompy.addToStudyInFather(solid_rev, grp_w, "wall_ch")
    # Injector fluid wall group
grp_w_of = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_w_of, faces_w_o+faces_w_f)
geompy.addToStudyInFather(solid_rev, grp_w_of, "wall_of")
    # Injector remaining wall group
grp_w_i = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_w_i, faces_w_i)
geompy.addToStudyInFather(solid_rev, grp_w_i, "wall_i")
    # Injector inlet oxidizer
grp_in_o = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_in_o, faces_in_o)
geompy.addToStudyInFather(solid_rev, grp_in_o, "inlet_o")
    # Injector inlet fuel
grp_in_f = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_in_f, faces_in_f)
geompy.addToStudyInFather(solid_rev, grp_in_f, "inlet_f")
    # Chamber outlet group
grp_out = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_out, faces_out)
geompy.addToStudyInFather(solid_rev, grp_out, "outlet")
    # 0° periodic face group
grp_per_0 = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_per_0, [faces_p_0])
geompy.addToStudyInFather(solid_rev, grp_per_0, "periodic_0")
    # x° periodic face group
grp_per_rev = geompy.CreateGroup(solid_rev, geompy.ShapeType["FACE"])
geompy.UnionList(grp_per_rev, [faces_p_rev])
geompy.addToStudyInFather(solid_rev, grp_per_rev, "periodic_rev")

#%% Setup mesh

def set_mesh(type, max_size, min_size, growth_rate, group = None):
    mesh_handle = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D, geom=group) if type == "2D" else mesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    params_handle = mesh_handle.Parameters()
    params_handle.SetFineness(smeshBuilder.Custom)
    params_handle.SetMaxSize(max_size)
    params_handle.SetMinSize(min_size)
    params_handle.SetGrowthRate(growth_rate)
    params_handle.SetUseSurfaceCurvature(False)
    params_handle.SetSecondOrder(False)
    return mesh_handle

# Initialize mesh
mesh = smeshpy.Mesh(solid_rev)

min_fine_size = (rp * lc) / 5.0
max_coarse_size = cf*lc

# Global element settings
mesh_3D = set_mesh("3D", max_coarse_size, min_fine_size, 0.2)
# 2D surface mesh refinement near injector
mesh_2D_inj = set_mesh("2D", rp*lc, min_fine_size, 0.15, grp_w_of)
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
    faces_id = [geompy.GetSubShapeID(solid_rev, face) for face in faces]
    return faces_id

# Chamber wall inflation layers
visc_w = add_viscous_layers(4.03e-4, 15, 1.01, get_faces_id(faces_w))

# Injector fluid wall inflation layers
visc_w_of = add_viscous_layers(2.17e-5, 15, 1.13, get_faces_id(faces_w_f+faces_w_o))

# Injector fluid wall inflation layers
visc_w_i = add_viscous_layers(4.03e-4, 15, 1.01, get_faces_id(faces_w_i))

#%% Compute mesh and export

is_done = mesh.Compute()

# Dictionary mapping OpenFOAM patch names to geompy groups
patch_groups = {
    "wall_ch": grp_w,
    "wall_of": grp_w_of,
    "wall_i": grp_w_i,
    "inlet_o": grp_in_o,
    "inlet_f": grp_in_f,
    "outlet": grp_out,
    "periodic_0": grp_per_0,
    "periodic_rev": grp_per_rev
}

for patch_name, geom_grp in patch_groups.items():
    mesh.GroupOnGeom(geom_grp, patch_name, SMESH.FACE)

if is_done and mesh_export:
    output_path = str(Path.cwd() / "run_openfoam" / "constant" / "polyMesh")
    # Export mesh directly
    exportToFoam(mesh, output_path)
    print(f"✅ Exported OpenFOAM Mesh {len(patch_groups)} patches to {output_path}")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()