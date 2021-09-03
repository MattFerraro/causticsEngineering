module CausticsEngineering

using DocStringExtensions

using Meshes, Images, Plots
gr()

include("parameters.jl")
include("utilities.jl")
include("plots.jl")
include("io.jl")

include("create_mesh.jl")



"""
$(SIGNATURES)
"""
function main()
    @assert size(ARGS) == (1,) "Intented usage is: julia create_mesh.jl image.png"

    img = Images.load(ARGS[1])
    return engineer_caustics(img)
end


"""
$(SIGNATURES)
"""
main()


export FaceMesh, main, engineer_caustics

end
