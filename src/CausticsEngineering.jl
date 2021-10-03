module CausticsEngineering

using DocStringExtensions
using Random, Images

using Plots
gr()

using Meshes, FileIO, MeshIO


include("../examples/personal/original.jl")


include("parameters.jl")
include("utilities.jl")
include("math_utilities.jl")
include("mesh.jl")
include("io.jl")
include("plots.jl")

include("create_mesh.jl")



"""
$(SIGNATURES)
"""
function main()
    @assert size(ARGS) == (1,) "Intented usage is: julia create_mesh.jl image.png"

    img = Images.load(ARGS[1])
    return engineer_caustics(img)
end

export
    # Constants
    Meters_Per_Pixel,
    Bottom_Offset,
    Top_Offset,

    # Types
    Vertex3D,
    FieldVertex3D,
    FaceMesh,

    # math
    average,
    average_absolute,
    laplacian,
    ∇,

    # Mesh
    create_solid,
    create_solid_as_dict,
    get_lens_pixels_area,

    # Plotting
    plot_as_quiver,

    # I/O
    save_obj!,

    # Main procedures
    propagate_poisson!,
    march_mesh!,
    engineer_caustics,
    main


end # End module
