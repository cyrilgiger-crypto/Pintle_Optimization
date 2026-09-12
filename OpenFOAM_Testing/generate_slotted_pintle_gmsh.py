#%% Imports
import gmsh
import numpy as np

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
cf = 3          # [-], coarsening factor outlet

#%% Define points

tand = lambda deg: np.tan(np.deg2rad(deg))
# sind = lambda deg: np.sin(np.deg2rad(deg))

a_bl = (np.pi/2 - np.deg2rad(th_b))/2   # blend half angle

gmsh.initialize()
gmsh.model.add(name)

points1 = []
points1.append(gmsh.model.geo.addPoint(L_po+L_op+t_pt, 0, 0, lc))
points1.append(gmsh.model.geo.addPoint(L_do, 0, 0, cf*lc))
points1.append(gmsh.model.geo.addPoint(L_do, D_do/2, 0, cf*lc))
points1.append(gmsh.model.geo.addPoint(0, D_do/2, 0, cf*lc))
points1.append(gmsh.model.geo.addPoint(0, D_po/2+t_an, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(-L_ex, D_po/2+t_an, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(-L_ex, D_po/2, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(L_po, D_po/2, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(L_po, D_po/2-t_po, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(L_po - (D_po/2-t_po-D_mp/2)/tand(th_t), D_mp/2, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(-L_ex, D_mp/2, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(-L_ex, D_pr/2, 0, rp*lc))
points1.append(gmsh.model.geo.addPoint(L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2, 0, rp*lc))

point_bl = gmsh.model.geo.addPoint(L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl), D_pr/2+r_bl, 0, rp*lc)

points2 = []
points2.append(gmsh.model.geo.addPoint(L_po+L_op-(D_po-D_pr)/2*tand(th_b)-r_bl*np.tan(a_bl)+r_bl*np.sin(2*a_bl),D_pr/2+r_bl*(1-np.cos(2*a_bl)), 0, rp*lc))
points2.append(gmsh.model.geo.addPoint(L_po+L_op, D_po/2, 0, rp*lc))
points2.append(gmsh.model.geo.addPoint(L_po+L_op+t_pt, D_po/2, 0, rp*lc))
points2.append(points1[0])

#%% Helper functions for inflation layer construction

def add_2d_inflation_layer(curve_tags, hwall_n, ratio, thickness):
    fan_points = set()
    for c in curve_tags:
        bnd = gmsh.model.getBoundary([(1, c)])
        for dim, ptag in bnd:
            fan_points.add(ptag)
    bl_field = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl_field, "CurvesList", curve_tags)
    gmsh.model.mesh.field.setNumber(bl_field, "hwall_n", hwall_n)
    gmsh.model.mesh.field.setNumber(bl_field, "ratio", ratio)
    gmsh.model.mesh.field.setNumber(bl_field, "thickness", thickness)
    gmsh.model.mesh.field.setNumbers(bl_field, "FanPointsList", list(fan_points))
    gmsh.model.mesh.field.setAsBoundaryLayer(bl_field)
    return bl_field

#%% Create domain

# create sub-curves
l1 = [gmsh.model.geo.addLine(points1[i], points1[i+1]) for i in range(len(points1) - 1)]
l2 = [gmsh.model.geo.addLine(points2[i], points2[i+1]) for i in range(len(points2) - 1)]
f1 = [gmsh.model.geo.addCircleArc(points1[-1], point_bl, points2[0])]

# connect all curves to create loop
c_all = l1 + f1 + l2
boundary = gmsh.model.geo.addCurveLoop(c_all)

# consturct surface from loop
surface = gmsh.model.geo.addPlaneSurface([boundary])
gmsh.model.geo.synchronize()

# apply 2D inflation layers to revolve
l_wall = [3,4,15]; l_inj_o = np.r_[8:11, 12:15, 16].tolist(); l_inj_f = [5, 7]
# bl1 = add_2d_inflation_layer(l_wall, 1e-5, 1.2, 1e-3)
# gmsh.model.mesh.field.setNumbers(bl1, "FanPointsList", [1, 2, 3, 5])
# bl2 = add_2d_inflation_layer(l_inj_f+l_inj_o, 1e-6, 1.2, 1e-4)
# gmsh.model.mesh.field.setNumbers(bl2, "FanPointsList", [5, 6, 7, 11, 12, 17])
# gmsh.model.mesh.setRecombine(2, surface)

# revolve
rev_out = gmsh.model.geo.revolve([(2,surface)], 0, 0, 0, 1, 0, 0, np.pi/2, recombine=True)
vol = [tag for dim, tag in rev_out if dim == 3][0]

# build geometry
gmsh.model.geo.synchronize()

#%% Tag surfaces

# find created surface tags from revolve
surf_rev = [tag for dim, tag in rev_out if dim==2]

# dict to map line tag to surface tag
lin2sur = dict(zip(c_all[1:], surf_rev[1:]))    # exclude line1 since its the center line

gmsh.model.addPhysicalGroup(2, [lin2sur[6]], name="inlet_fl")
gmsh.model.addPhysicalGroup(2, [lin2sur[11]], name="inlet_ox")
gmsh.model.addPhysicalGroup(2, [lin2sur[2]], name="outlet")
gmsh.model.addPhysicalGroup(1, [1], name="axis")

tag_wall = [lin2sur[i] for i in l_wall]
gmsh.model.addPhysicalGroup(2, tag_wall, name="wall")

tag_wall_f = [lin2sur[i] for i in l_inj_f]
tag_wall_o = [lin2sur[i] for i in l_inj_o]
gmsh.model.addPhysicalGroup(2, tag_wall_f, name="wall_f")
gmsh.model.addPhysicalGroup(2, tag_wall_o, name="wall_o")

gmsh.model.addPhysicalGroup(2, [surface], name="periodic_0")
gmsh.model.addPhysicalGroup(2, [surf_rev[0]], name="periodic_rev")

gmsh.model.addPhysicalGroup(3, [vol], name="fluid")

#%% Mesh geometry and run gmsh

gmsh.model.mesh.generate(2)
gmsh.model.mesh.generate(3)

gmsh.fltk.run()

gmsh.finalize()

#%%