using Revise, Debugger

using Images
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

using CausticsEngineering
mesh, imageBW = engineer_caustics(image);


# Check a few values to make sure they make sense
mesh.corners.r
minimum(mesh.corners.vr)
maximum(mesh.corners.vr)

mesh.corners.c
minimum(mesh.corners.vc)
maximum(mesh.corners.vc)

minimum(mesh.corners.ϕ)
maximum(mesh.corners.ϕ)
mean_ϕ = sum(mesh.corners.ϕ) / length(mesh.corners.ϕ)

# Plot the last vector field
plot_as_quiver(
    mesh,
    stride = 20,
    scale = 1.0,
    max_length = 20,
    flipxy = true,
    reverser = false,
    reversec = false,
)


# Generate dictionaries and arrays containing the vertices and triangles.
triangle_dict, vertex_dict, triangle_index, vertex_index =
    create_solid(mesh; bottom_distance = Bottom_Offset, top_distance = Top_Offset)


# Can be done in 2 steps
#
# triangle_dict, vertex_dict = create_solid_as_dict(
#     mesh;
#     bottom_distance = Bottom_Offset,
#     top_distance = Top_Offset
# )
#
# # Generate array or coordinates and indices usable to create traditional mesh objects
# _, _, triangle_index, vertex_index = create_solid(
#     triangle_dict,
#     vertex_dict;
#     bottom_distance = Bottom_Offset,
#     top_distance = Top_Offset,
#     )


# Save as an obj file
save_obj!(
    triangle_index,
    vertex_index,
    "./examples/result_cat_mesh.obj",
    scale = Float64(Meters_Per_Pixel),
    scaleh = Float64(Meters_Per_Pixel),
)


# Convert to a mesh structure
using Meshes

# 2D Mesh
vertex_index2D = [(v[1], v[2]) for v ∈ vertex_index]
vertices2 = Meshes.Point2.(vertex_index2D)
cat_mesh2D = Meshes.SimpleMesh(vertices, connections)

# 3D Mesh
vertices = Meshes.Point3.(vertex_index)
connections = Meshes.connect.(triangle_index)
cat_mesh = Meshes.SimpleMesh(vertices, connections)


# NOT WORKING...
using MeshViz

# Explicitly import the preferred backend (depending on use case)
import CairoMakie
MeshViz.viz(cat_mesh)
MeshViz.viz(cat_mesh2D)
