using Revise, Debugger

using Images
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

using CausticsEngineering
mesh, image_caustics = engineer_caustics(image);
