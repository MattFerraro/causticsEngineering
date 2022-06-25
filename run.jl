using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Images, CausticsEngineering

const input_file = "./data/input.jpg"

image = Images.load(input_file) # Check current working directory with pwd()
engineer_caustics(image);
