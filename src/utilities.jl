"""
$(TYPEDEF)

## Coordinates

"""
struct Vertex3D
    x::Float64
    y::Float64
    z::Float64

    vx::Float64
    vy::Float64

    Vertex3D(x, y, z, vx, vy) = new(x, y, z, vx, vy)
    Vertex3D() = new(0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
$(SIGNATURES)
"""
function dist(p1::Tuple{Float64,Float64,Float64}, p2::Tuple{Float64,Float64,Float64})
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    dz = p2[3] - p1[3]
    return sqrt(dx^2 + dy^2 + dz^2)
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
function area(
    p1::Tuple{Float64,Float64,Float64},
    p2::Tuple{Float64,Float64,Float64},
    p3::Tuple{Float64,Float64,Float64},
)
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a + b + c) / 2.0

    return sqrt(s * (s - a) * (s - b) * (s - c))
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

Centroid of three points.
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
function fill_borders!(mat::AbstractMatrix{T}, val::T) where {T}
    mat[begin:end, begin] .= val
    mat[begin:end, end] .= val
    mat[begin, begin:end] .= val
    mat[end, begin:end] .= val
end


"""
$(TYPEDEF)

## Coordinates

Origin is top left, going right and down .`x` goes horizontal and follows columns. x is within 1 and width
`y` goes vertical and follows rows. y is within 1 and height.

# Velocity

Velocity vector along `x` and `y` (velocity along z not necessary).

Implementation is a struc of arrays. Easier to vectorise.

"""
mutable struct FieldVertex3D
    size::Tuple{Int,Int}

    x::AbstractMatrix{Float64}
    y::AbstractMatrix{Float64}
    z::AbstractMatrix{Float64}

    vx::AbstractMatrix{Float64}
    vy::AbstractMatrix{Float64}

    function FieldVertex3D(height, width)
        # mx = zeros(Float64, height + 1, width + 1)
        # my = zeros(Float64, height + 1, width + 1)
        # mz = zeros(Float64, height + 1, width + 1)

        # mvx = zeros(Float64, height + 1, width + 1)
        # mvy = zeros(Float64, height + 1, width + 1)

        mx = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        my = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        mz = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000

        mvx = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        mvy = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000

        new((height, width), mx, my, mz, mvx, mvy)
    end
end

Base.size(fv::FieldVertex3D) = FieldVertex3D.size


"""
$(SIGNATURES)
"""
function Vertex3D(fv::FieldVertex3D, row, col)
    height, width = size(fv)

    if (1 <= row <= height + 1) && (1 <= col <= width + 1)
        Vertex3D(
            fv.x[row, col],
            fv.y[row, col],
            fv.z[row, col],
            fv.vx[row, col],
            fv.vy[row, col],
        )
    else
        missing
    end
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
struct FaceMesh
    size::Tuple{Int,Int}
    topleft::FieldVertex3D

    toptriangles::AbstractMatrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}
    bottriangles::AbstractMatrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}

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

    m_topleft = FieldVertex3D(height, width)
    t_triangles =
        Matrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}(undef, height, width)
    b_triangles =
        Matrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}(undef, height, width)

    for row ∈ 1:height, col ∈ 1:width
        t_triangles[row, col] = ((row, col), (row + 1, col), (row, col + 1))
        b_triangles[row, col] = ((row, col + 1), (row + 1, col), (row + 1, col + 1))
    end

    return FaceMesh((height, width), m_topleft, t_triangles, b_triangles)
end


"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
Base.size(mesh::FaceMesh) = mesh.size


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

Converts a triangle as a triplet of references to mesh vertices to a triplet of 3D coordinates.
"""
function top_triangle3D(mesh::FaceMesh, row::Int, col::Int)
    height, width = size(mesh)
    if 1 <= row <= height + 1 && 1 <= col <= width + 1
        t = mesh.toptriangles[row, col]
        t1 = CartesianIndex(t[1][1], t[1][2])
        t2 = CartesianIndex(t[2][1], t[2][2])
        t3 = CartesianIndex(t[3][1], t[3][2])
        p1 = Vertex3D(
            mesh.topleft.x[t1],
            mesh.topleft.y[t1],
            mesh.topleft.z[t1],
            mesh.topleft.vx[t1],
            mesh.topleft.vy[t1],
        )
        p2 = Vertex3D(
            mesh.topleft.x[t2],
            mesh.topleft.y[t2],
            mesh.topleft.z[t2],
            mesh.topleft.vx[t2],
            mesh.topleft.vy[t2],
        )
        p3 = Vertex3D(
            mesh.topleft.x[t3],
            mesh.topleft.y[t3],
            mesh.topleft.z[t3],
            mesh.topleft.vx[t3],
            mesh.topleft.vy[t3],
        )
        return (p1, p2, p3)
    else
        return missing
    end
end

top_triangle3D(mesh::FaceMesh, ci::CartesianIndex{2}) = top_triangle3D(mesh, ci[1], ci[2])


function bot_triangle3D(mesh::FaceMesh, row::Int, col::Int)
    height, width = size(mesh)
    if 1 <= row <= height + 1 && 1 <= col <= width + 1
        t = mesh.bottriangles[row, col]
        t1 = CartesianIndex(t[1][1], t[1][2])
        t2 = CartesianIndex(t[2][1], t[2][2])
        t3 = CartesianIndex(t[3][1], t[3][2])
        p1 = Vertex3D(
            mesh.topleft.x[t1],
            mesh.topleft.y[t1],
            mesh.topleft.z[t1],
            mesh.topleft.vx[t1],
            mesh.topleft.vy[t1],
        )
        p2 = Vertex3D(
            mesh.topleft.x[t2],
            mesh.topleft.y[t2],
            mesh.topleft.z[t2],
            mesh.topleft.vx[t2],
            mesh.topleft.vy[t2],
        )
        p3 = Vertex3D(
            mesh.topleft.x[t3],
            mesh.topleft.y[t3],
            mesh.topleft.z[t3],
            mesh.topleft.vx[t3],
            mesh.topleft.vy[t3],
        )
        return (p1, p2, p3)
    else
        return missing
    end
end

bot_triangle3D(mesh::FaceMesh, ci::CartesianIndex{2}) = bot_triangle3D(mesh, ci[1], ci[2])


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

Centroid of a specific mesh triangle.
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

Centroid of a specific mesh triangle.
"""
function centroid(mesh::FaceMesh, row::Int, col::Int; side = Union{:top,:bottom})
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
function area(mesh::FaceMesh, ci::CartesianIndex; side = Union{:top,:bottom})
    if side == :top
        return area(top_triangle3D(mesh, ci)...)
    elseif side == :bottom
        return area(top_triangle3D(mesh, ci)...)
    else
        return missing
    end
end

"""
$(SIGNATURES)
"""
function area(mesh::FaceMesh, row::Int, col::Int; side = Union{:top,:bottom})
    if side == :top
        return area(top_triangle3D(mesh, row, col)...)
    elseif side == :bottom
        return area(top_triangle3D(mesh, row, col)...)
    else
        return missing
    end
end
