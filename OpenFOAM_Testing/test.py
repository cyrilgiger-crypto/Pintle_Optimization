import gmsh
import math

gmsh.initialize()
gmsh.model.add("pintle_3d_sector")

# ---------------------------------------------------------
# Parametric Inputs
# ---------------------------------------------------------
n_holes = 12
sector_deg = 360.0 / n_holes
sector_rad = math.radians(sector_deg)

r_pintle = 0.005       # Pintle radius
r_outer = 0.1        # Outer orifice radius
l_chamber = 0.1      # Downstream chamber length
l_pintle_tip = 0.015   # Length of pintle body inside domain

lc = 0.0010            # Base mesh size

# ---------------------------------------------------------
# 1. Create 2D Profile in X-Y Plane
# ---------------------------------------------------------
p1 = gmsh.model.occ.addPoint(0.0, r_pintle, 0.0, lc)
p2 = gmsh.model.occ.addPoint(0.0, r_outer, 0.0, lc*5)
p3 = gmsh.model.occ.addPoint(l_chamber, r_outer, 0.0, lc*5)
p4 = gmsh.model.occ.addPoint(l_chamber, 0.0, 0.0, lc*5)
p5 = gmsh.model.occ.addPoint(l_pintle_tip, 0.0, 0.0, lc)
p6 = gmsh.model.occ.addPoint(l_pintle_tip, r_pintle, 0.0, lc)

l1 = gmsh.model.occ.addLine(p1, p2)
l2 = gmsh.model.occ.addLine(p2, p3)
l3 = gmsh.model.occ.addLine(p3, p4)
l4 = gmsh.model.occ.addLine(p4, p5)
l5 = gmsh.model.occ.addLine(p5, p6)
l6 = gmsh.model.occ.addLine(p6, p1)

wire = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4, l5, l6])
face = gmsh.model.occ.addPlaneSurface([wire])

# ---------------------------------------------------------
# 2. Revolve 2D Face Around X-Axis to create 3D Sector
# ---------------------------------------------------------
# revolve(dimTags, x, y, z, ax, ay, az, angle)
v = gmsh.model.occ.revolve([(2, face)], 0, 0, 0, 1, 0, 0, sector_rad)

gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Generate DD triangular mesh
gmsh.model.mesh.generate(3)
gmsh.write("pintle_3d.msh")

gmsh.fltk.run()

gmsh.finalize()

print(f"--> Revolved {sector_deg} deg 3D sector successfully!")