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
