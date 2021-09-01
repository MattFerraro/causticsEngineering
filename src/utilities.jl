"""
$(TYPEDEF)

## Coordinates

Origin is top left, going right and down .`x` goes horizontal and follows columns. x is within 1 and width
`y` goes vertical and follows rows. y is within 1 and height.

# Velocity

Velocity vector along `x` and `y` (velocity along z not necessary).

"""
mutable struct Vertex3D
    x::Float64
    y::Float64
    z::Float64

    vx::Float64
    vy::Float64

    Vertex3D() = new(0.0, 0.0, Height_Offset, 0.0, 0.0)
end


"""
$(SIGNATURES)
"""
function dist(p1::Vertex3D, p2::Vertex3D)
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dz = p2.z - p1.z
    return sqrt(dx^2 + dy^2 + dz^2)
end


"""
$(SIGNATURES)

A midpoint is the average between two points
"""
midpoint(p1::Vertex3D, p2::Vertex3D) =
    Vertex3D((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0, (p1.z + p2.z) / 2.0, 0.0, 0.0)


"""
$(SIGNATURES)

Barycentre of three points.
"""
centroid(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D) = Vertex3D(
    (p1.x + p2.x + p3.x) / 3.0,
    (p1.y + p2.y + p3.y) / 3.0,
    (p1.z + p2.z + p3.z) / 3.0,
    0,
    0,
)


"""
$(SIGNATURES)
"""
function triangle_area(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D)
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a + b + c) / 2.0

    return sqrt(s * (s - a) * (s - b) * (s - c))
end


"""
$(TYPEDEF)

The algorithm works by allocating a rectangle from the exit mesh (the face where the light exits the block
of acrylate). Initially, each rectangle representsa single pixel on that face and faces s single pixel on
the caustics picture.

Each rectangle is split into 2 trangles: one upper-left, the other bottom-right like:

```
*------*
|    / |
|   /  |
|  /   |
| /    |
*------*
```


Although the corners of each triangle match corners of the rectangles, it is sometimes easier to think in terms
of triangles rather than corners coming from different rectangles.
"""
Triangle = Tuple{UInt64,UInt64,UInt64}


"""
$(TYPEDEF)

A Rectagle mesh represents a single face of the block. The most important face is the one facing the cautics.
First the surface is split into rectangles each associated with a pixel. Then, each rectangle is broken
into a top-left triangle and a bottom-right triangle.

`VertexList` is the list of all the vertices (i.e. points) in the mesh. the list of triangles is split in 2.
Each triangle is a triplet of indices, each pointing to the relevant position in the list of vertices.

The vertices are numbered by rows, then by columns.
"""
struct RectangleMesh
    height::UInt64
    width::UInt64
    vertexList::Vector{Vertex3D}
    topTriangles::Matrix{Triangle}
    botTriangles::Matrix{Triangle}

    as_index::Function

    function RectangleMesh(height::Int, width::Int)

        to_index(r, c) = r + (c - 1) * (width + 1)

        vertices = Vector{Vertex3D}(undef, (height + 1) * (width + 1))
        for row = 1:height+1, col = 1:width+1
            vertices[to_index(row, col)] =
                Vertex3D(Float64(row), Float64(col), Height_Offset, 0.0, 0.0)
        end

        top_triangles = Matrix{Triangle}(undef, height, width)
        bot_triangles = Matrix{Triangle}(undef, height, width)
        for row = 1:height, col = 1:width
            top_triangles[row, col] =
                (to_index(row, col), to_index(row + 1, col), to_index(row, col + 1))
            bot_triangles[row, col] =
                (to_index(row, col + 1), to_index(row + 1, col), to_index(row + 1, col + 1))
        end

        return new(height, width, vertexList, top_triangles, bot_triangles, to_index)
    end
end



"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
function base.size(mesh::RectangleMesh)
    return mesh.height, mesh.width
end

"""
$(SIGNATURES)
"""
triangle3D(mesh::Mesh, t::Triangle) =
    [mesh.vertexList(t[1]), mesh.vertexList(t[1]), mesh.vertexList(t[1])]

"""
$(SIGNATURES)
"""
centroid(mesh::Mesh, t::Triangle) = centroid(triangle3D(mesh, t)...)


"""
$(SIGNATURES)
"""
triangle_area(mesh::Mesh, t::Triangle) = triangle_area(triangle3D(mesh, t)...)



function top_triangle(mesh::RectangleMesh, row::Int, col::Int)
    height, width = size(mesh)

    if (1 <= row <= height) && (1 <= col <= width)
        return mesh.topTriangles(row, col)
    else
        return nothing
    end
end

function bottom_triangle(mesh::RectangleMesh, row::Int, col::Int)
    height, width = size(mesh)

    if (1 <= row <= height) && (1 <= col <= width)
        return mesh.botTriangles(row, col)
    else
        return nothing
    end
end
