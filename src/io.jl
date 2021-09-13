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

"""
function save_obj!(
    triangle_index,
    vertex_index::AbstractMatrix{Float64},
    filename::String;
    scale = 1.0,
    scaleh = 1.0,
    reverse = false,
    flipxy = false,
)

    open(filename, "w") do io
        println(io, "solid engineered_caustics")

        xs = vertex_index[:, 1]
        ys = vertex_index[:, 2]
        zs = vertex_index[:, 3]
        for i ∈ 1:length(xs)
            if flipxy
                println(io, "v $(ys[i] * scale) $(xs[i][1] * scale) $(zs[i] * scaleh)")
            else
                println(io, "v $(xs[i] * scale) $(ys[i][1] * scale) $(zs[i] * scaleh)")
            end
        end

        t1 = triangle_index[1]
        t2 = triangle_index[2]
        t3 = triangle_index[3]
        for i ∈ 1:length(t1)
            println(io, "f $(t1[i]) $(t2[i]) $(t3[i])")
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
