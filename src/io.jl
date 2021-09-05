using FileIO, MeshIO


"""
$(SIGNATURES)

"""
function convert(::Meshes.SimpleMesh, mesh::FaceMesh)
    height, width = size(mesh)

    points = Meshes.Point3.([mesh.topleft[ci] for ci ∈ CartesianIndices(mesh.topleft)])

    top_connections =
        [connect(mesh.toptriangles[ci]) for ci ∈ CartesianIndices(mesh.topleft)]
    bot_connections =
        [connect(mesh.bottriangles[ci]) for ci ∈ CartesianIndices(mesh.topleft)]

    return SimpleMesh(points, vcat(top_connections, bot_connections))
end



"""
$(SIGNATURES)

This function saves the mesh object in stl format.

The format difinition is sourced from [https://en.wikipedia.org/wiki/STL_(file_format)]().

TO REFACTOR.
"""
function save_stl!(
    mesh::FaceMesh,
    filename::String;
    scale = 1.0,
    scalez = 1.0,
    reverse = false,
    flipxy = false,
)

    return

    height, width = size(mesh)

    open(filename, "w") do io
        println(io, "solid engineered_caustics")

        for row ∈ 1:height, col ∈ 1:width
            # Top triangle
            top_triangle = top_triangle(mesh, row, col)
            n = centroid(top_triangle)

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            # Bottom triangle
            bot_triangle = bottom_triangle(mesh, row, col)
            n = centroid(bot_triangle)

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end
        end

        # CHECK what dims exactly represents. Number of triangles?
        println(io, "dims $(2*mesh.width) $(2*mesh.height)")

        println(io, "endsolid engineered_caustics")
    end
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function stl_2_mesh(filename)
    lines = readlines(filename)

    vertexLines = [l for l in lines if startswith(l, "v")]
    nodeList = Vector{Vertex3D}(undef, size(vertexLines))

    count = 0
    for line in vertexLines
        count += 1
        elements = split(line, " ")
        x = parse(Float64, elements[2])
        y = parse(Float64, elements[3])
        z = parse(Float64, elements[4]) * 10
        pt = Vertex3D(x, y, z, 0, 0)
        nodeList[count] = pt
    end

    faceLines = [l for l in lines if startswith(l, "f")]
    triangles = Vector{Triangle}(undef, size(faceLines))
    for line in faceLines
        elements = split(line, " ")
        triangles[line] = Triangle(
            parse(Int64, elements[2]),
            parse(Int64, elements[3]),
            parse(Int64, elements[4]),
        )
    end

    dimsLines = [l for l in lines if startswith(l, "dims")]
    elements = split(dimsLines[1], " ")

    return RectangleMesh(
        nodeList,
        triangles,
        parse(Int64, elements[2]),
        parse(Int64, elements[3]),
    )
end
