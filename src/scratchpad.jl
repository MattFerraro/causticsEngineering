using Revise, Debugger

using Images
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

using CausticsEngineering
mesh, imageBW = engineer_caustics(image);

CausticsEngineering.plot_as_quiver(
    mesh,
    stride = 20,
    scale = 1.0,
    max_length = 20,
    flipxy = true,
    reverser = false,
    reversec = false,
)



triangle_dict, vertex_dict =
    create_solid_as_dict(mesh; bottom_distance = Bottom_Offset, top_distance = Top_Offset)

_, _, triangle_index, vertex_index = create_solid(
    triangle_dict,
    vertex_dict;
    bottom_distance = Bottom_Offset,
    top_distance = Top_Offset,
)


CausticsEngineering.save_stl!(
    triangle_index,
    vertex_index,
    "./examples/result_mesh.obj",
    scale = Float64(Meters_Per_Pixel),
    scaleh = Float64(Meters_Per_Pixel),
)
