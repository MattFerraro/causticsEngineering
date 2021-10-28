if size(ARGS) == (0,)
    println("Intented usage is: julia run.jl image.png")
    file = "examples/cat_posing.jpg"
else
    file = ARGS[1]
end
println("Working on file $(file)")

using Pkg
Pkg.activate(".")

using Images, CausticsEngineering

img = load(file)
return engineer_caustics(img, string(file[1:end-3], "obj"))
