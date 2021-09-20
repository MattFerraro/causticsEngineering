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
image = Images.load("./examples/personal/portrait.jpg"); # Check current working directory with pwd()
image = Images.load("./examples/personal/bilal.jpg"); # Check current working directory with pwd()

mesh, imageBW = engineer_caustics(image);

mesh, imageBW = CausticsEngineering.original_engineer_caustics(image);


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









imageBW = Float64.(Gray.(image))

# Set global size parameters
height, width = size(imageBW)
global N_Pixel_Height = height
global N_Pixel_Width = width
global Meters_Per_Pixel = Caustics_Long_Side / N_Pixel_Height
println("Image size: $(N_Pixel_Height) x $(N_Pixel_Width)")

# mesh is the same size as the image with an extra row/column to have coordinates to
# cover each image corner with a triangle.
mesh = FaceMesh(N_Pixel_Height, N_Pixel_Width)

# The total energy going through the lens is equal to the amount of energy on the caustics
total_energy_lens = N_Pixel_Height * N_Pixel_Width * 1.0     # 1 unit of energy per pixel
total_energy_caustics = sum(imageBW)
average_energy_per_pixel = average(imageBW)

# imageBW is normalised to the same (sort of) _energy_ as the original image.
imageBW = imageBW / average_energy_per_pixel

lens_pixels_area = get_lens_pixels_area(mesh)

# Positive error => the triangle needs to shrink (less light)
error_luminosity = Float64.(lens_pixels_area - imageBW)
error_luminosity = error_luminosity .- average(error_luminosity)

# Save the loss image as a png
println(
    """
Luminosity:
    Max/min pixel areas: $(maximum(lens_pixels_area)) / $(minimum(lens_pixels_area))
    Average pixel areas: $(average_absolute(lens_pixels_area))
    Average pixel areas: $(average(lens_pixels_area))

    Max/min luminosity error: $(maximum(error_luminosity)) / $(minimum(error_luminosity))
    Average abs error: $(average_absolute(error_luminosity))
    Average error: $(average(error_luminosity))
    """,
)



##########################################################
## STEPPING

using Revise, Debugger, Images, Plots;
gr();

using CausticsEngineering


# engineer_caustics
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()
imageBW = Float64.(Gray.(image));
imageBW /= average(imageBW);
sum(imageBW)

height, width = size(imageBW)




# solve_velocity_potential
mesh = CausticsEngineering.FaceMesh(height, width);
CausticsEngineering.field_summary(mesh.corners.ϕ)

origmesh = CausticsEngineering.squareMesh(width + 1, height + 1);
origϕ = zeros(width, height);
CausticsEngineering.field_summary(origϕ)


# height, width = size(mesh.corners.ϕ)
area_distorted_corners = CausticsEngineering.get_lens_pixels_area(mesh);
origarea_distorted_corners = CausticsEngineering.getPixelArea(origmesh);
CausticsEngineering.field_summary(area_distorted_corners)
CausticsEngineering.field_summary(area_distorted_corners - origarea_distorted_corners)

error_luminosity = Float64.(area_distorted_corners - imageBW);
CausticsEngineering.field_summary(error_luminosity)

origerror_luminosity = Float64.(origarea_distorted_corners - imageBW);
CausticsEngineering.field_summary(origerror_luminosity)

CausticsEngineering.field_summary(error_luminosity - origerror_luminosity)


error_luminosity = error_luminosity .- average(error_luminosity);
CausticsEngineering.field_summary(error_luminosity)

origerror_luminosity = origerror_luminosity .- average(origerror_luminosity);
CausticsEngineering.field_summary(origerror_luminosity)

CausticsEngineering.field_summary(error_luminosity - origerror_luminosity)



mesh = CausticsEngineering.FaceMesh(height, width);
CausticsEngineering.field_summary(mesh.corners.ϕ)
Lϕ = CausticsEngineering.laplacian(mesh.corners.ϕ);
CausticsEngineering.field_summary(Lϕ)
δ = Lϕ - error_luminosity;
δ .*= CausticsEngineering.ω / 4.0;
CausticsEngineering.field_summary(δ)


mesh.corners.ϕ[1:end-1, 1:end-1] .+= δ;
mesh.corners.ϕ .-= average(mesh.corners.ϕ);
CausticsEngineering.field_summary(mesh.corners.ϕ)


mesh = CausticsEngineering.FaceMesh(height, width);
CausticsEngineering.field_summary(CausticsEngineering.laplacian(mesh.corners.ϕ) - error_luminosity)

max_update, Lϕ, δ = CausticsEngineering.propagate_poisson!(mesh.corners.ϕ, error_luminosity);
CausticsEngineering.field_summary(mesh.corners.ϕ)
CausticsEngineering.field_summary(Lϕ)
CausticsEngineering.field_summary(δ)


origϕ = zeros(width, height);
max_update, origLϕ, origδ = CausticsEngineering.orig_propagate_poisson!(origϕ, origerror_luminosity);
CausticsEngineering.field_summary(origϕ)
CausticsEngineering.field_summary(origLϕ)
CausticsEngineering.field_summary(origδ)







mesh = CausticsEngineering.FaceMesh(height, width);
new_update = 10_000
old_update = 2*new_update
for i ∈ 1:1_000
    old_update = new_update
    new_update, Lϕ, δ = CausticsEngineering.propagate_poisson!(mesh.corners.ϕ, error_luminosity)
    i % 50 == 0 && println("""
        $(i) => $(max_update)
            $(CausticsEngineering.field_summary(mesh.corners.ϕ))
            $(CausticsEngineering.field_summary(Lϕ))
            $(CausticsEngineering.field_summary(δ))
            """
    )
    if new_update > old_update || i == 5_000
        new_update, Lϕ, δ = CausticsEngineering.orig_propagate_poisson!(origϕ, origerror_luminosity)
        println("""
        $(i) => $(new_update)
            $(CausticsEngineering.field_summary(mesh.corners.ϕ))
            $(CausticsEngineering.field_summary(Lϕ))
            $(CausticsEngineering.field_summary(δ))
            """
            )
        break
    end
end

origϕ = zeros(width, height);
for i ∈ 1:5_000
    max_update,origLϕ, origδ = CausticsEngineering.orig_propagate_poisson!(origϕ, origerror_luminosity);
    i % 500 == 0 && println("""
        $(i) => $(max_update)
            $(CausticsEngineering.field_summary(origϕ))
            $(CausticsEngineering.field_summary(origLϕ))
            $(CausticsEngineering.field_summary(origδ))
            """
    )
    if max_update <= 1e-3
        max_update, origLϕ, origδ = CausticsEngineering.orig_propagate_poisson!(origϕ, origerror_luminosity);
        println(CausticsEngineering.field_summary(origϕ))
        println(CausticsEngineering.field_summary(origLϕ))
        println(CausticsEngineering.field_summary(origδ))
        break
    end
end










CausticsEngineering.save_plot_scalar_field!(error_luminosity, "trace_1", imageBW)

max_update = Inf
old_update = Inf

any(isnan.(error_luminosity))

any(isnan.(mesh.corners.r))
any(isnan.(mesh.corners.c))
any(isnan.(mesh.corners.vr))
any(isnan.(mesh.corners.vc))

any(isnan.(mesh.corners.ϕ))

i = 0
begin
    i += 1
    old_update = max_update
    max_update = CausticsEngineering.propagate_poisson!(mesh.corners.ϕ, error_luminosity)
    # i % 100 == 0 && println(max_update)
    # (abs(old_update - max_update) <= 1e-6 || abs(max_update) <= 1e-4) && break
end

max_update
old_update


f()


∇ϕᵤ, ∇ϕᵥ = CausticsEngineering.∇(mesh.corners.ϕ)

sum(∇ϕᵤ)
average(∇ϕᵤ)
average_absolute(∇ϕᵤ)


height, width = size(imageBW)
mesh.corners.vr .= -∇ϕᵤ
mesh.corners.vc .= -∇ϕᵥ
list_triangles = vcat(
    [
        CausticsEngineering.triangle3D(mesh, row, col, :top) for row ∈ 1:height, col ∈ 1:width
    ],
    [
        CausticsEngineering.triangle3D(mesh, row, col, :bottom) for row ∈ 1:height, col ∈ 1:width
    ],
)

list_maximum_t = [time for time ∈ CausticsEngineering.find_maximum_t.(list_triangles) if !isnothing(time) && time > 0.]
minimum(list_maximum_t)
maximum(list_maximum_t)

δ = minimum(list_maximum_t) / 2.0

mesh.corners.r += δ * ∇ϕᵤ
mesh.corners.c += δ * ∇ϕᵥ


end







nothing;
