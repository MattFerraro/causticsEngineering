using Images, CausticsEngineering


image = Images.load("./examples/cat_posing.jpg") # Check current working directory with pwd()
mesh, image_caustics = engineer_caustics(image);
nothing;

image_caustics'
