module CausticsEngineering

using DocStringExtensions

using Images, Plots
gr()

include("parameters.jl")
include("utilities.jl")
include("create_mesh.jl")

export main, engineer_caustics

end
