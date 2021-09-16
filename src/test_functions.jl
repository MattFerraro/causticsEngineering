using CausticsEngineering
using Test

"""
$(SIGNATURES)

"""
function testSquareMesh!()
    mesh = create_mesh(100, 50)

    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodeArray[1, 1].x = 8
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodes[1].y += 12
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function testSolidify!()
    println("Testing solidification")
    width = 100
    height = 100
    origMesh = create_mesh(width, height)

    for y = 1:height, x = 1:width
        x2 = (x - width / 2) / width
        y2 = (y - height / 2) / height
        value = x2 * x2 + y2 * y2
        origMesh.nodeArray[x, y].z = 15 - value * 25
    end

    save_stl!(origMesh, "./examples/testSolidify.obj")
    solidMesh = solidify(origMesh, 0)
    save_stl!(solidMesh, "./examples/testSolidify2.obj")
end





function testAreaTriangle(x1, y1, x2, y2, x3, y3)
    a = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    b = sqrt((x3 - x2)^2 + (y3 - y2)^2)
    c = sqrt((x1 - x3)^2 + (y1 - y3)^2)
    s = (a + b + c) / 2.0

    surface_sq = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16.0
    return surface_sq <= 1e-100 ? 0.0 : surface_sq
end

function doTestAreaTriangle(N = 10_000)
    X1 = 10_000 .* rand(N)
    Y1 = 10_000 .* rand(N)
    X2 = 10_000 .* rand(N)
    Y2 = 10_000 .* rand(N)
    X3 = 10_000 .* rand(N)
    Y3 = 10_000 .* rand(N)

    for x1 ∈ X1, y1 ∈ Y1, x2 ∈ X2, y2 ∈ Y2, x3 ∈ X3, y3 ∈ Y3

        _ = testAreaTriangle(x1, y1, x2, y2, x3, y3) > 0.0
    end
end

@testset "triangle area" begin
    for x1 ∈ 10_000 .* rand(10_000),
        y1 ∈ 10_000 .* rand(10_000),
        x2 ∈ 10_000 .* rand(10_000),
        y2 ∈ 10_000 .* rand(10_000),
        x3 ∈ 10_000 .* rand(10_000),
        y3 ∈ 10_000 .* rand(10_000)

        @test testAreaTriangle(x1, y1, x2, y2, x3, y3) >= 0.0
    end
end

@testset "More triangle area" begin
    @test area(0, 1, 1 , 0, 0, 0)
end


v1 = CausticsEngineering.Vertex3D(0.0, 1.0, 0.0 , 0.0, 0.0)
v2 = CausticsEngineering.Vertex3D(1.0, 0.0, 0.0 , 0.0, 0.0)
v3 = CausticsEngineering.Vertex3D(0.0, 0.0, 0.0 , 0.0, 0.0)
CausticsEngineering.area(v1, v2, v3)

v1 = CausticsEngineering.Vertex3D(0.0, 2.0, 0.0 , 0.0, 0.0)
v2 = CausticsEngineering.Vertex3D(2.0, 0.0, 0.0 , 0.0, 0.0)
v3 = CausticsEngineering.Vertex3D(0.0, 0.0, 0.0 , 0.0, 0.0)
CausticsEngineering.area(v1, v2, v3)

fm = CausticsEngineering.FaceMesh(512, 512);
fm.toptriangles[1, 1]
CausticsEngineering.triangle3D(fm, 1, 1, :top)
CausticsEngineering.area(CausticsEngineering.triangle3D(fm, 1, 1, :top)...)
CausticsEngineering.area(CausticsEngineering.triangle3D(fm, 1, 1, :bottom)...)

CausticsEngineering.area(CausticsEngineering.triangle3D(fm, 126, 13, :top)...)
CausticsEngineering.area(CausticsEngineering.triangle3D(fm, 126, 13, :bottom)...)



v1 = [7250.768146661309, 4849.994680605491]
v2 = [1725.2665442814496, 8169.783790398901]
v3 = [6560.1725812433, 5264.912964907599]

a = sqrt(sum((v2 - v1) .^ 2))
b = sqrt(sum((v3 - v2) .^ 2))
c = sqrt(sum((v1 - v3) .^ 2))
s = (a + b + c) / 2.0

s * (s - a) * (s - b) * (s - c)
(a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) / 16.0

b + c
a

10^8 * (b + c - a)
