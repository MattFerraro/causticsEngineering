"""
$(TYPEDEF)

## Coordinates

"""
struct Vertex3D
    r::Float64
    c::Float64
    z::Float64

    vr::Float64
    vc::Float64

    Vertex3D(r, c, z, vr, vc) = new(r, c, z, vr, vc)
    Vertex3D() = new(0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
$(SIGNATURES)
"""
function dist(p1::Tuple{Float64,Float64,Float64}, p2::Tuple{Float64,Float64,Float64})
    dr = p2[1] - p1[1]
    dc = p2[2] - p1[2]
    dz = p2[3] - p1[3]
    return sqrt(dr^2 + dc^2 + dz^2)
end

"""
$(SIGNATURES)
"""
function dist(p1::Vertex3D, p2::Vertex3D)
    dr = p2.r - p1.r
    dc = p2.c - p1.c
    dz = p2.z - p1.z
    return sqrt(dr^2 + dc^2 + dz^2)
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
    Vertex3D((p1.r + p2.r) / 2.0, (p1.c + p2.c) / 2.0, (p1.z + p2.z) / 2.0, 0.0, 0.0)


"""
$(SIGNATURES)

Centroid of three points.
"""
centroid(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D) = Vertex3D(
    (p1.r + p2.r + p3.r) / 3.0,
    (p1.c + p2.c + p3.c) / 3.0,
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

    r::AbstractMatrix{Float64}
    y::AbstractMatrix{Float64}
    z::AbstractMatrix{Float64}

    vr::AbstractMatrix{Float64}
    vc::AbstractMatrix{Float64}

    function FieldVertex3D(height, width)
        # mr = zeros(Float64, height + 1, width + 1)
        # my = zeros(Float64, height + 1, width + 1)
        # mz = zeros(Float64, height + 1, width + 1)

        # mvr = zeros(Float64, height + 1, width + 1)
        # mvc = zeros(Float64, height + 1, width + 1)

        rows = repeat(Float64.(1:height+1), 1, width + 1)
        cols = repeat(Float64.(1:width+1)', height + 1, 1)

        mr = rows + rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        my = cols + rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        mz = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000

        mvr = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
        mvc = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000

        new((height, width), mr, my, mz, mvr, mvc)
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
            fv.r[row, col],
            fv.y[row, col],
            fv.z[row, col],
            fv.vr[row, col],
            fv.vc[row, col],
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

    FaceMesh(s, tl, tt, bt) = new(s, tl, tt, bt)

    function FaceMesh(height::Int, width::Int)
        m_topleft = FieldVertex3D(height, width)
        t_triangles = Matrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}(
            undef,
            height,
            width,
        )
        b_triangles = Matrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}(
            undef,
            height,
            width,
        )

        for row ∈ 1:height, col ∈ 1:width
            t_triangles[row, col] = ((row, col), (row + 1, col), (row, col + 1))
            b_triangles[row, col] = ((row, col + 1), (row + 1, col), (row + 1, col + 1))
        end

        return FaceMesh((height, width), m_topleft, t_triangles, b_triangles)
    end
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
    if 1 <= row <= height && 1 <= col <= width
        t = mesh.toptriangles[row, col]
        t1_row, t1_col = t[1]
        t2_row, t2_col = t[2]
        t3_row, t3_col = t[3]

        p1 = Vertex3D(
            mesh.topleft.r[t1_row, t1_col],
            mesh.topleft.c[t1_row, t1_col],
            mesh.topleft.z[t1_row, t1_col],
            mesh.topleft.vr[t1_row, t1_col],
            mesh.topleft.vc[t1_row, t1_col],
        )
        p2 = Vertex3D(
            mesh.topleft.r[t2_row, t2_col],
            mesh.topleft.c[t2_row, t2_col],
            mesh.topleft.z[t2_row, t2_col],
            mesh.topleft.vr[t2_row, t2_col],
            mesh.topleft.vc[t2_row, t2_col],
        )
        p3 = Vertex3D(
            mesh.topleft.r[t3_row, t3_col],
            mesh.topleft.c[t3_row, t3_col],
            mesh.topleft.z[t3_row, t3_col],
            mesh.topleft.vr[t3_row, t3_col],
            mesh.topleft.vc[t3_row, t3_col],
        )
        return (p1, p2, p3)
    else
        return missing
    end
end

top_triangle3D(mesh::FaceMesh, ci::CartesianIndex{2}) = top_triangle3D(mesh, ci[1], ci[2])


function bot_triangle3D(mesh::FaceMesh, row::Int, col::Int)
    height, width = size(mesh)
    if 1 <= row <= height && 1 <= col <= width
        t = mesh.bottriangles[row, col]
        t1_row, t1_col = t[1]
        t2_row, t2_col = t[2]
        t3_row, t3_col = t[3]

        p1 = Vertex3D(
            mesh.topleft.r[t1_row, t1_col],
            mesh.topleft.c[t1_row, t1_col],
            mesh.topleft.z[t1_row, t1_col],
            mesh.topleft.vr[t1_row, t1_col],
            mesh.topleft.vc[t1_row, t1_col],
        )
        p2 = Vertex3D(
            mesh.topleft.r[t2_row, t2_col],
            mesh.topleft.c[t2_row, t2_col],
            mesh.topleft.z[t2_row, t2_col],
            mesh.topleft.vr[t2_row, t2_col],
            mesh.topleft.vc[t2_row, t2_col],
        )
        p3 = Vertex3D(
            mesh.topleft.r[t3_row, t3_col],
            mesh.topleft.c[t3_row, t3_col],
            mesh.topleft.z[t3_row, t3_col],
            mesh.topleft.vr[t3_row, t3_col],
            mesh.topleft.vc[t3_row, t3_col],
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
