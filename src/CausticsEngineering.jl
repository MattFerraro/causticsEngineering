module CausticsEngineering

using DocStringExtensions
using Random, Meshes, Images

using Plots
gr()

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

    # Mesh
    create_solid,
    create_solid_as_dict,

    # Plotting
    plot_as_quiver,

    # I/O
    save_obj!,

    # Main procedures
    main,
    engineer_caustics


end # End module
