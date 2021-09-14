using Revise, Debugger

using Images, Plots;
gr();
using CausticsEngineering

image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

image = Images.load("./examples/personal/salvador_dali_1.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/salvador_dali_2.jpg"); # Check current working directory with pwd()

image = Images.load("./examples/personal/statue_of_liberty_1.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/statue_of_liberty_2.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/statue_of_liberty_3.jpg"); # Check current working directory with pwd()

image = Images.load("./examples/personal/caricature.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/image.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/bilal.jpg"); # Check current working directory with pwd()

mesh, imageBW = engineer_caustics(image; clamp_correction = true);


imageBW = Float64.(Gray.(image));
average_energy_caustics = CausticsEngineering.average(imageBW)
imageBW = imageBW / average_energy_caustics;

(sum(imageBW) / length(imageBW) - 1.0) * 1e12


# Check a few values to make sure they make sense
mesh.corners.r
mesh.corners.vr
minimum(mesh.corners.vr)
maximum(mesh.corners.vr)

mesh.corners.c
mesh.corners.vc
minimum(mesh.corners.vc)
maximum(mesh.corners.vc)

mesh.corners.ϕ
minimum(mesh.corners.ϕ)
maximum(mesh.corners.ϕ)
mean_ϕ = sum(mesh.corners.ϕ) / length(mesh.corners.ϕ)

# Plot the last vector field
p = plot_as_quiver(mesh, n_steps = 60, scale = 5.0, max_length = 20)



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

# Area of all the distorted pixels
am = CausticsEngineering.get_area_corners(mesh)
maximum(am)
minimum(am)

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

# Zone to plot
row_min = 100.0
row_max = 200.0
col_min = 100.0
col_max = 200.0

function clip_trangle_index(
    triangle_index,
    vertex_index;
    row_min = 100.0,
    row_max = 200.0,
    col_min = 100.0,
    col_max = 200.0,
)

    function is_in_zone(t)::Bool
        p1, p2, p3 = t
        v1_r, v1_c = vertex_index[p1]
        v2_r, v2_c = vertex_index[p2]
        v3_r, v3_c = vertex_index[p3]

        return row_min <= v1_r <= row_max &&
               col_min <= v1_c <= col_max &&
               row_min <= v2_r <= row_max &&
               col_min <= v2_c <= col_max &&
               row_min <= v3_r <= row_max &&
               col_min <= v3_c <= col_max
    end

    return [t for t ∈ triangle_index if is_in_zone(t)]
end

# Select relevant triangles
triangle_index3D = clip_trangle_index(triangle_index, vertex_index)
connections = Meshes.connect.(triangle_index3D)

# 3D Mesh
vertices3D = Meshes.Point3.(vertex_index)
cat_mesh3D = Meshes.SimpleMesh(vertices3D, connections)

# 2D Mesh
vertex_index2D = [(v[1], v[2]) for v ∈ vertex_index]
vertices2D = Meshes.Point2.(vertex_index2D)
cat_mesh2D = Meshes.SimpleMesh(vertices2D, connections)


# NOT WORKING...
using MeshViz

# Explicitly import the preferred backend (depending on use case)
import GLMakie
MeshViz.viz(cat_mesh3D)
MeshViz.viz(cat_mesh2D)
