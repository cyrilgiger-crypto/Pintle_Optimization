import salome
salome.salome_init()

import GEOM
from salome.geom import geomBuilder # type: ignore
import SMESH
from salome.smesh import smeshBuilder # type: ignore

# --- 1. GEOMETRY CREATION ---
geompy = geomBuilder.New()

dx, dy, dz = 100.0, 100.0, 50.0
box = geompy.MakeBoxDXDYDZ(dx, dy, dz)

# Identify bottom (z = 0) and top (z = dz) faces
bottom_face = geompy.GetFaceNearPoint(box, geompy.MakeVertex(dx / 2.0, dy / 2.0, 0.0))
top_face = geompy.GetFaceNearPoint(box, geompy.MakeVertex(dx / 2.0, dy / 2.0, dz))

bottom_face_id = geompy.GetSubShapeID(box, bottom_face)

geompy.addToStudy(box, 'Cuboid')
geompy.addToStudyInFather(box, bottom_face, 'Bottom_Face')
geompy.addToStudyInFather(box, top_face, 'Top_Face')

# --- 2. GLOBAL MESH SETUP ---
smesh = smeshBuilder.New()
mesh = smesh.Mesh(box)

# Global 3D Tetrahedral Mesher
netgen = mesh.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
netgen_params = netgen.Parameters()
netgen_params.SetMaxSize(10.0)
netgen_params.SetMinSize(2.0)

# --- 3. LOCAL SUB-MESHES (FACE REFINEMENT) ---

# Sub-mesh 1: Bottom Face -> 5 mm max size
submesh_bottom = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D, geom=bottom_face)
params_bottom = submesh_bottom.Parameters()
params_bottom.SetMaxSize(5.0)
params_bottom.SetMinSize(1.0)

# Sub-mesh 2: Top Face -> 10 mm max size
submesh_top = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D, geom=top_face)
params_top = submesh_top.Parameters()
params_top.SetMaxSize(10.0)
params_top.SetMinSize(2.0)

# --- 4. VISCOUS INFLATION LAYER SETUP ---
viscous_layers = smesh.CreateHypothesis('ViscousLayers')
viscous_layers.SetTotalThickness(5.0)
viscous_layers.SetNumberLayers(5)
viscous_layers.SetStretchFactor(1.2)
viscous_layers.SetFaces([bottom_face_id], 0)  # Extrude inward on bottom face

mesh.AddHypothesis(viscous_layers)

# --- 5. COMPUTE MESH ---
is_done = mesh.Compute()

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()