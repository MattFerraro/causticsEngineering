module CausticsEngineering

using DocStringExtensions

using meshes

using Images, Plots
gr()

include("parameters.jl")
include("utilities.jl")
include("create_mesh.jl")

export main, engineer_caustics

end
