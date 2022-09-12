# Work in a temporary environment
using Pkg
Pkg.activate(; temp = true)

# Speed up by avoiding updating the repository when adding packages
Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true

# Add useful package
Pkg.add([
    "Revise", "Images"
])

Pkg.develop(path = @__DIR__)

using Revise, Images

using CausticsEngineering

# Check current working directory with pwd()
image = Images.load("./examples/personal/caricature.jpg")
engineer_caustics(image);
