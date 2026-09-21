#%% Imports
import math
import os
import numpy as np
from pathlib import Path
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
L_po = 3.0e-3   # [m], pintle post length
r_fl = 1.0e-3   # [m], fillet radius for fuel exit for smooth BL growth and prevent Co number spikes

# Oxidizer side
D_pr = 3.0e-3   # [m], pintle rod diameter
D_mp = 6.0e-3   # [m], minimum internal diamater of pintle post
th_t = 20.0     # [°], top interal pintle angle
th_b = 20.0     # [°], bottom interal pintle angle
L_op = 0.5e-3   # [m], pintle slit opening
t_po = 0.5e-3   # [m], post thickness at slit
t_pt = 1e-3     # [m], pintle plate thickness
r_bl = 2e-3     # [m], blend radius on rod

# Domain dimensions
L_ex = 10e-3    # [m], extension of ox and f channels in -x
D_do = 70e-3    # [m], domain external wall diameter
L_do = 60e-3   # [m], domain length

# Mesh settings
lc = 1.0e-3     # [m], base mesh size
rp = 0.2        # [-], refinement factor near pintle
cf = 3.0        # [-], coarsening far away
angle = 90      # [°], revolution angle for periodic BC

# Misc
mesh_export = False

#%% Define geometry points

tand = lambda deg: np.tan(np.deg2rad(deg))
a_bl = (np.pi/2 - np.deg2rad(th_b))/2   # blend half angle

points1 = [
    (L_po+L_op+t_pt, 0, 0),
    (L_do, 0, 0),
    (L_do, D_do/2, 0),
    (0, D_do/2, 0),
    (0, D_po/2+t_an+r_fl, 0)]

point_fl_c = (-r_fl, D_po/2+t_an+r_fl, 0)

points2 = [
    (-r_fl, D_po/2+t_an, 0),
    (-L_ex, D_po/2+t_an, 0),
    (-L_ex, D_po/2, 0),
    (L_po, D_po/2, 0),
    (L_po, D_po/2-t_po, 0),
    (L_po - (D_po/2-t_po-D_mp/2)/tand(th_t), D_mp/2, 0),
    (-L_ex, D_mp/2, 0),
    (-L_ex, D_pr/2, 0),
    (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2, 0)]

point_bl_c = (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2+r_bl, 0)

points3 = [
    (L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl)+r_bl*np.sin(2*a_bl),D_pr/2+r_bl*(1-np.cos(2*a_bl)), 0),
    (L_po+L_op, D_po/2, 0),
    (L_po+L_op+t_pt, D_po/2, 0)]

#%% Build geometry

# Define vertices
vert1 = [geompy.MakeVertex(*pt) for pt in points1]
v_center1 = geompy.MakeVertex(*point_fl_c)
vert2 = [geompy.MakeVertex(*pt) for pt in points2]
v_center2 = geompy.MakeVertex(*point_bl_c)
vert3 = [geompy.MakeVertex(*pt) for pt in points3]

# Connect first part with lines
lin1 = [geompy.MakeLineTwoPnt(vert1[i], vert1[i+1]) for i in range(len(vert1) - 1)]
arc1 = [geompy.MakeArcCenter(v_center1, vert1[-1], vert2[0], False)]
lin2 = [geompy.MakeLineTwoPnt(vert2[i], vert2[i+1]) for i in range(len(vert2) - 1)]
arc2  = [geompy.MakeArcCenter(v_center2, vert2[-1], vert3[0], False)]
lin3 = [geompy.MakeLineTwoPnt(vert3[i], vert3[i+1]) for i in range(len(vert3) - 1)]
closing_line = [geompy.MakeLineTwoPnt(vert3[-1], vert1[0])]

# Make closed surface
wire = geompy.MakeWire(lin1+arc1+lin2+arc2+lin3+closing_line)
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
lin_w = [lin1[i] for i in [2, 3]]
lin_w_f = arc1 + [lin2[i] for i in [0, 2]]
lin_w_o = [lin2[i] for i in [3, 4, 5, 7]] + arc2 + lin3
lin_w_i = closing_line
    # In and outlet lines
lin_in_f = [lin2[1]]
lin_in_o = [lin2[6]]
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

#%% Export stl and add patch naming

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

# Set export path
output_dir = os.getcwd() / Path("run_openfoam_hex_amr/constant/triSurface")
output_dir.mkdir(parents=True, exist_ok=True)
combined_stl = output_dir / f"{name}.stl"

deflection = 1.0e-5

with open(combined_stl, "w") as outfile:
    for patch_name, grp in patch_groups.items():
        temp_stl = output_dir / f"temp_{patch_name}.stl"

        geompy.ExportSTL(grp, str(temp_stl), True, deflection, True)

        with open(temp_stl, "r") as infile:
            for line in infile:
                if line.strip().startswith("solid"):
                    outfile.write(f"solid {patch_name}\n")
                elif line.strip().startswith("endsolid"):
                    outfile.write(f"endsolid {patch_name}\n")
                else:
                    outfile.write(line)

        temp_stl.unlink()

print(f"Multi-patch STL created successfully at: {combined_stl}")