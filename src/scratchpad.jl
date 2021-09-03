using Revise

using Images
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

using CausticsEngineering
mesh, image_caustics = engineer_caustics(image);

Gray.(image_caustics')


typeof(CartesianIndex(3, 5)[1])

t = Matrix{Tuple{CartesianIndex(2),CartesianIndex(2),CartesianIndex(2)}}(undef, 2, 2)

t[1, 1] = (CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2))
typeof((CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2)))

f = FaceMesh(5, 5);

f.toptriangles[CartesianIndex(5, 5)]

f.topleft


v3 = CausticsEngineering.top_triangle3D(f, CartesianIndex(5, 5))
CausticsEngineering.area(v3...)

v3 = CausticsEngineering.bot_triangle3D(f, CartesianIndex(5, 5))
CausticsEngineering.area(v3...)

CausticsEngineering.get_pixel_area(f)

for _ = 1:100
    h = rand(Float64, (10, 10))
    d = rand(Float64, (10, 10))
    print(CausticsEngineering.relax(h, d))
end

c = t[1]

f.topleft(c)

t = CausticsEngineering.top_triangle3D(f, CartesianIndex(1, 1))
