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
"""
function dist(p1::Vec3f, p2::Vec3f)
    v = p2 .- p1
    return sqrt(sum(v .^ 2))
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
Rectangle = Nothing


"""
$(TYPEDEF)
"""
struct FaceMesh
    size::Tuple{UInt,UInt}
    topleft::AbstractMatrix{Vertex3D}
    toptriangles::AbstractMatrix{Tuple{UInt,UInt,UInt}}
    bottriangles::AbstractMatrix{Tuple{UInt,UInt,UInt}}

    as_index::Function

    function FaceMesh(height::Int, width::Int)

        function as_index(r, c)
            if 1 <= r <= height + 1 && 1 <= c <= width
                return r + (c - 1) * (width + 1)
            else
                return nothing
            end
        end

        m_topleft = Matrix{Vertex3D}(undef, height, width)
        t_triangles = Matrix{Tuple{UInt,UInt,UInt}}(undef, height, width)
        b_triangles = Matrix{Tuple{UInt,UInt,UInt}}(undef, height, width)
        for ci ∈ CartesianIndices(m_topleft)
            row = ci[1]
            col = ci[2]

            m_topleft .= (row, col, Top_Offset, 0.0, 0.0)

            t_triangles[ci] =
                (as_index(row, col), as_index(row + 1, col), as_index(row, col + 1))
            b_triangles[ci] =
                (as_index(row, col + 1), as_index(row + 1, col), as_index(row + 1, col + 1))
        end

        return new((height, width), m_topleft, t_triangles, b_triangles)
    end
end

function convert(Meshes.SimpleMesh, mesh::FaceMesh)
    height, width = size(mesh)

    points = Meshes.Point3.([mesh.topleft[ci] for ci ∈ CartesianIndices(mesh.topleft)])

    top_connections =
        [connect(mesh.toptriangles[ci]) for ci ∈ CartesianIndices(mesh.topleft)]
    bot_connections =
        [connect(mesh.bottriangles[ci]) for ci ∈ CartesianIndices(mesh.topleft)]

    return SimpleMesh(points, vcat(top_connections, bot_connections))
end


"""
$(TYPEDEF)
"""
struct Triangle
    mesh::FaceMesh
    points::Tuple{UInt,UInt,UInt}

    function Triangle(mesh::FaceMesh)
        return Triangle(
            mesh,
            Tuple(mesh.as_index(1, 1), mesh.as_index(2, 1), mesh.as_index(1, 2)),
        )
    end
end


"""
$(SIGNATURES)
"""
function area(p1::Vec3f, p2::Vec3f, p3::Vec3f)
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a + b + c) / 2.0

    return sqrt(s * (s - a) * (s - b) * (s - c))
end


"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
base.size(mesh::FaceMesh) = mesh.size


"""
$(SIGNATURES)

Converts a triangle as a triplet of references to mesh vertices to a triplet of 3D coordinates.
"""
function top_triangle3D(mesh::Mesh, ci::CartesianIndex)
    if ci ∈ CartesianIndices(mesh.topleft)
        t = mesh.top_triangle[ci]
        p1 = mesh.topleft[t[1]]
        p2 = mesh.topleft[t[2]]
        p3 = mesh.topleft[t[3]]
        return (p1, p2, p3)
    else
        return missing
    end
end

top_triangle3D(mesh::Mesh, row::UInt, width::UInt) =
    top_triangle3D(mesh, CartesianIndex(row, col))

function bot_triangle3D(mesh::Mesh, ci::CartesianIndex)
    if ci ∈ CartesianIndices(mesh.topleft)
        t = mesh.bot_triangle[ci]
        p1 = mesh.topleft[t[1]]
        p2 = mesh.topleft[t[2]]
        p3 = mesh.topleft[t[3]]
        return (p1, p2, p3)
    else
        return missing
    end
end

bot_triangle3D(mesh::Mesh, row::UInt, width::UInt) =
    bot_triangle3D(mesh, CartesianIndex(row, col))

function triangle3D(mesh::Mesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return top_triangle3D(mesh, ci)
    elseif side == :bottom
        return bot_triangle3D(mesh, ci)
    else
        return missing
    end
end

triangle3D(mesh::Mesh, row::Uint, col::UInt; side = Union{:top,:bottom}) =
    triangle3D(mesh, CartesianIndex(row, col); side)


"""
$(SIGNATURES)
"""
function centroid(mesh::Mesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return centroid(top_triangle3D(mesh, row, col)...)
    elseif side == :bottom
        return centroid(bot_triangle3D(mesh, row, col)...)
    else
        return missing
    end
end


"""
$(SIGNATURES)
"""
function area(mesh::Mesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return area(top_triangle3D(mesh, ci)...)
    elseif side == :bottom
        return area(top_triangle3D(mesh, ci)...)
    else
        return missing
    end
end
