using StaticArrays

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

    Vertex3D(x, y, z, vx, vy) = new(x, y, z, vx, vy)
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
function area(v1::Vertex3D, v2::Vertex3D, v3::Vertex3D)
    a = dist(v1, v2)
    b = dist(v2, v3)
    c = dist(v3, v1)
    s = (a + b + c) / 2.0

    return sqrt(s * (s - a) * (s - b) * (s - c))
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
    size::Tuple{Int,Int}
    topleft::AbstractMatrix{Vertex3D}
    toptriangles::AbstractMatrix{Tuple{CartesianIndex,CartesianIndex,CartesianIndex}}
    bottriangles::AbstractMatrix{Tuple{CartesianIndex,CartesianIndex,CartesianIndex}}

    as_index::Function

    FaceMesh(s, tl, tt, bt) = new(s, tl, tt, bt)
end

function FaceMesh(height::Int, width::Int)

    function as_index(r, c)
        # Adding "+ 1" since it is the size of the `topleft` matrix
        if 1 <= r <= height + 1 && 1 <= c <= width + 1
            return r + (c - 1) * (width + 1)
        else
            return nothing
        end
    end
    as_index(ci::CartesianIndex) = as_index(ci[1], ci[2])

    m_topleft = Matrix{Vertex3D}(undef, height + 1, width + 1)
    t_triangles = Matrix{Tuple{CartesianIndex{2},CartesianIndex{2},CartesianIndex{2}}}(
        undef,
        height,
        width,
    )
    b_triangles = Matrix{Tuple{CartesianIndex{2},CartesianIndex{2},CartesianIndex{2}}}(
        undef,
        height,
        width,
    )


    for ci ∈ CartesianIndices(m_topleft)
        row = ci[1]
        col = ci[2]
        m_topleft[ci] = Vertex3D(row, col, Top_Offset, 0.0, 0.0)
    end

    for ci ∈ CartesianIndices(t_triangles)
        row = ci[1]
        col = ci[2]

        t_triangles[ci] = (
            CartesianIndex(row, col),
            CartesianIndex(row + 1, col),
            CartesianIndex(row, col + 1),
        )
        b_triangles[ci] = (
            CartesianIndex(row, col + 1),
            CartesianIndex(row + 1, col),
            CartesianIndex(row + 1, col + 1),
        )
    end

    return FaceMesh((height, width), m_topleft, t_triangles, b_triangles)
end


"""
$(TYPEDEF)
"""
struct Triangle
    mesh::FaceMesh
    points::Tuple{Int,Int,Int}

    function Triangle(mesh::FaceMesh)
        return Triangle(
            mesh,
            Tuple(mesh.as_index(1, 1), mesh.as_index(2, 1), mesh.as_index(1, 2)),
        )
    end
end



"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
Base.size(mesh::FaceMesh) = mesh.size


"""
$(SIGNATURES)

Converts a triangle as a triplet of references to mesh vertices to a triplet of 3D coordinates.
"""
function top_triangle3D(mesh::FaceMesh, ci::CartesianIndex)
    if ci ∈ CartesianIndices(mesh.toptriangles)
        t = mesh.toptriangles[ci]
        p1 = mesh.topleft[t[1]]
        p2 = mesh.topleft[t[2]]
        p3 = mesh.topleft[t[3]]
        return (p1, p2, p3)
    else
        return missing
    end
end

top_triangle3D(mesh::FaceMesh, row::Int, width::Int) =
    top_triangle3D(mesh, CartesianIndex(row, col))

function bot_triangle3D(mesh::FaceMesh, ci::CartesianIndex)
    if ci ∈ CartesianIndices(mesh.bottriangles)
        t = mesh.bottriangles[ci]
        p1 = mesh.topleft[t[1]]
        p2 = mesh.topleft[t[2]]
        p3 = mesh.topleft[t[3]]
        return (p1, p2, p3)
    else
        return missing
    end
end

bot_triangle3D(mesh::FaceMesh, row::Int, width::Int) =
    bot_triangle3D(mesh, CartesianIndex(row, col))


function triangle3D(mesh::FaceMesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return top_triangle3D(mesh, ci)
    elseif side == :bottom
        return bot_triangle3D(mesh, ci)
    else
        return missing
    end
end

triangle3D(mesh::FaceMesh, row::Int, col::Int; side = Union{:top,:bottom}) =
    triangle3D(mesh, CartesianIndex(row, col); side)


"""
$(SIGNATURES)
"""
function centroid(mesh::FaceMesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return centroid(top_triangle3D(mesh, ci)...)
    elseif side == :bottom
        return centroid(bot_triangle3D(mesh, ci)...)
    else
        return missing
    end
end


"""
$(SIGNATURES)
"""
function area(mesh::FaceMesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return area(top_triangle3D(mesh, ci)...)
    elseif side == :bottom
        return area(top_triangle3D(mesh, ci)...)
    else
        return missing
    end
end
