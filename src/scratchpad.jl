using Revise, Debugger

using Images
image = Images.load("./examples/cat_posing.jpg"); # Check current working directory with pwd()

using CausticsEngineering
mesh, imageBW = engineer_caustics(image);

triangle_dict, vertex_dict, triangle_index, vertex_index = create_solid(mesh)

save_stl!(
    triangle_index, vertex_index,
    "./examples/result_mesh.obj",
    scale = Float64(Meters_Per_Pixel),
    scaleh = Float64(Meters_Per_Pixel),
)









Gray.(imageBW)

mesh.pixel.r[1:10, 1:10]
mesh.pixel.c[1:10, 1:10]
mesh.pixel.ϕ[1:10, 1:10]
mesh.pixel.vr[1:10, 1:10]
mesh.pixel.vc[1:10, 1:10]

mesh.toptriangles[1:10]
mesh.bottriangles[1:10]

mesh.pixel.r[250:260, 250:260]
mesh.pixel.c[250:260, 250:260]
mesh.pixel.ϕ[250:260, 250:260]


using Test

t3 = (
    CausticsEngineering.Vertex3D(5.0, 5.0, 0.01, 0.1, -2.0),
    CausticsEngineering.Vertex3D(6.0, 5.0, 0.01, 0.1, 2.0),
    CausticsEngineering.Vertex3D(5.5, 5.5, 0.01, 0.5, -3.0),
)

CausticsEngineering.find_maximum_t(t3)


typeof(CartesianIndex(3, 5)[1])

t = Matrix{Tuple{CartesianIndex(2),CartesianIndex(2),CartesianIndex(2)}}(undef, 2, 2)

t[1, 1] = (CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2))
typeof((CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2)))


f.toptriangles[CartesianIndex(5, 5)]

f = FaceMesh(5, 5);
v3 = CausticsEngineering.top_triangle3D(f, CartesianIndex(5, 5))
CausticsEngineering.find_maximum_t(v3)

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

f.pixel(c)

t = CausticsEngineering.top_triangle3D(f, CartesianIndex(1, 1))



# Indexing structs
struct TT
    a::Int
    b::Int
end

using StructArrays


aa = StructArrays()

aa[:].a = 0


[1, 2, 3] |> filter(isnan)


begin
    r1, c1, ϕ1, vr1, vc1 = (rand(5) .- 0.5) .* 20
    r2, c2, ϕ2, vr2, vc2 = (rand(5) .- 0.5) .* 20
    r3, c3, ϕ3, vr3, vc3 = (rand(5) .- 0.5) .* 20

    vertex1 = Vertex3D(r1, c1, ϕ1, vr1, vc1)
    vertex2 = Vertex3D(r2, c2, ϕ2, vr2, vc2)
    vertex3 = Vertex3D(r3, c3, ϕ3, vr3, vc3)

    println(
        CausticsEngineering.area(vertex1, vertex2, vertex3),
        "   ",
        CausticsEngineering.find_maximum_t(vertex1, vertex2, vertex3),
    )
end


r1, c1, ϕ1, vr1, vc1 = (rand(5) .- 0.5) .* 20
r2, c2, ϕ2, vr2, vc2 = (rand(5) .- 0.5) .* 20
r3, c3, ϕ3, vr3, vc3 = (rand(5) .- 0.5) .* 20

p1 = Vertex3D(r1, c1, ϕ1, vr1, vc1)
p2 = Vertex3D(r2, c2, ϕ2, vr2, vc2)
p3 = Vertex3D(r3, c3, ϕ3, vr3, vc3)

# To make the calculation simpler, everything is translated so that A is at the
# origin of the plane and its velocity is nil.
Br = p2.r - p1.r
Bc = p2.c - p1.c
Cr = p3.r - p1.r
Cc = p3.c - p1.c

t_vBr = p2.vr - p1.vr
t_vBc = p2.vc - p1.vc
t_vCr = p3.vr - p1.vr
t_vCc = p3.vc - p1.vc

# After this, given that Ar = Ac = t_vAr = t_vAc = 0, the area is nil iff
# (Br + t_vBr) (Cc + t_vCc ) - (Cr + t_vCr) (Bc + t_vBc) = 0.
# After expansion and reshuffling to have a quadratic equation where t
# is the variable, the coefficients of that equation are:
a = t_vCc * t_vBr - t_vBc * t_vCr
b = -Bc * t_vCr - Cr * t_vBc + Br * t_vCc + Cc * t_vBr
c = Br * Cc - Cr * Bc

discriminant = b^2 - 4a * c
d = discriminant >= 0.0 ? sqrt(discriminant) : 0
t1 = (-b - d) / 2a
t2 = (-b + d) / 2a
t1, t2, CausticsEngineering.smallest_positive((-b - d) / 2a, (-b + d) / 2a)


Vertex3D()




a = Vertex3D(
    489.20336199266757,
    190.09282140997482,
    -58505.12856934825,
    -39628.63636211195,
    -50804.04679723241,
)
b = Vertex3D(
    490.30273062001396,
    189.42556748453276,
    29305.470378387537,
    -10020.45529953039,
    48181.962585623834,
)
c = Vertex3D(
    489.89918140614986,
    189.67049879639754,
    -18876.492207236297,
    39648.76478266787,
    -65706.41461026447,
)
CausticsEngineering.find_maximum_t(a, b, c)

a = Vertex3D(30.0, 271.0, 0.8354160335261099, -0.053188093456838725, -0.00538632666184502)
b = Vertex3D(31.0, 270.0, 0.8353725644270783, -0.005828514543554375, -0.05323156255587036)
c = Vertex3D(
    30.499999997149644,
    270.5,
    0.8886041269829487,
    0.04690021252675358,
    0.04690021225938901,
)
CausticsEngineering.find_maximum_t(a, b, c)

ab = CausticsEngineering.dist(a, b)
ac = CausticsEngineering.dist(a, c)
cb = CausticsEngineering.dist(c, b)

s = (ab + ac + cb) / 2.0

s - ab
s - ab
s - ab

[[0 1 0], [1 -4 1], [0 1 0]]


    flow_up = zeros(Float64, size(∇²ϕ))
    flow_down = zeros(Float64, size(∇²ϕ))
    flow_left = zeros(Float64, size(∇²ϕ))
    flow_right = zeros(Float64, size(∇²ϕ))

    flow_up[1:end, 1:end] .= padded_ϕ[1:end-2, 2:end-1] - padded_ϕ[2:end-1, 2:end-1]
    flow_down[1:end, 1:end] .= padded_ϕ[3:end, 2:end-1] - padded_ϕ[2:end-1, 2:end-1]
    flow_left[1:end, 1:end] .= padded_ϕ[2:end-1, 1:end-2] - padded_ϕ[2:end-1, 2:end-1]
    flow_right[1:end, 1:end] .= padded_ϕ[2:end-1, 3:end] - padded_ϕ[2:end-1, 2:end-1]

    fill_borders!(flow_up, 0.0)
    fill_borders!(flow_down, 0.0)
    fill_borders!(flow_left, 0.0)
    fill_borders!(flow_right, 0.0)




end

ϕ

k = ([[0.0 1.0 0.0]; [1.0 -4.0 1.0]; [0.0 1.0 0.0]])

m = [[0.0 1.0 0.0 4.0]; [1.0 -4.0 1.0 2.0]; [0.0 1.0 0.0 3.0]; [1.0 -4.0 1.0 2.0]]
p = zeros(Float64, 4, 4)


p[1, 1] = sum(m[1:3, 1:3] .* k)
p

vertex_dict = Dict{String, Tuple{Float64, Float64, Float64}}()

name = "a"
vertex_dict[name] = (1.0, 1.0, 1.0)
vertex_dict