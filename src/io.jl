using FileIO, MeshIO


# Utility function to name a vertex
vertex_name(face::String, row::Int64, col::Int64)::String = "Vertex_$(name)_$(row)_$(col)"

# Utility functions to create the names
triangle_name(face::String, row::Int64, col::Int64, side::Union{Symbol,Symbol})::String =
    "Tri_$(face)_$(row)_$(col)_$(side)"



"""
$(SIGNATURES)

"""
function convert(::Meshes.SimpleMesh, mesh::FaceMesh)
    height, width = size(mesh)

    points = Meshes.Point3.([mesh.corners[ci] for ci ∈ CartesianIndices(mesh.corners)])

    top_connections =
        [connect(mesh.toptriangles[ci]) for ci ∈ CartesianIndices(mesh.corners)]
    bot_connections =
        [connect(mesh.bottriangles[ci]) for ci ∈ CartesianIndices(mesh.corners)]

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
    scaleh = 1.0,
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
                println(io, "v $(n.c * scale) $(n.r * scale) $(n.h * scaleh)")
            else
                println(io, "v $(n.r * scale) $(n.c * scale) $(n.h * scaleh)")
            end

            # Bottom triangle
            bot_triangle = bottom_triangle(mesh, row, col)
            n = centroid(bot_triangle)

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

"""
function save_stl!(
    triangle_index::AbstractVector{},
    vertex_index::AbstractVector{},
    filename::String;
    scale = 1.0,
    scaleh = 1.0,
    reverse = false,
    flipxy = false,
)

    open(filename, "w") do io
        println(io, "solid engineered_caustics")

        for i ∈ 1:length(vertex_index)
            if flipxy
                println(
                    io,
                    "v $(vertex_index[2] * scale) $(vertex_index[1] * scale) $(vertex_index[3] * scaleh)",
                )
            else
                println(
                    io,
                    "v $(vertex_index[1] * scale) $(vertex_index[2] * scale) $(vertex_index[3] * scaleh)",
                )
            end
        end

        for i ∈ 1:length(triangle_index)
            println(
                io,
                "f $(triangle_index[i][1]) $(triangle_index[i][2]) $(triangle_index[i][3])",
            )
        end

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
        r = parse(Float64, elements[2])
        c = parse(Float64, elements[3])
        h = parse(Float64, elements[4]) * 10
        pt = Vertex3D(r, c, h, 0, 0)
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
