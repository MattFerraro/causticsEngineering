"""
$(TYPEDEF)

## Coordinates

"""
struct Vertex3D
    r::Float64
    c::Float64
    h::Float64

    vr::Float64
    vc::Float64

    Vertex3D(r, c, h, vr, vc) = new(r, c, h, vr, vc)
    Vertex3D() = new(0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
$(SIGNATURES)
"""
function dist(p1::Tuple{Float64,Float64,Float64}, p2::Tuple{Float64,Float64,Float64})
    dr = p2[1] - p1[1]
    dc = p2[2] - p1[2]
    dh = p2[3] - p1[3]
    return sqrt(dr^2 + dc^2 + dh^2)
end

"""
$(SIGNATURES)
"""
function dist(p1::Vertex3D, p2::Vertex3D)
    dr = p2.r - p1.r
    @assert 0.0 <= abs(dr) <= 1e6 "Distance (row coordinate = $(dr)) between $(p1) and $(p2) makes no sense "

    dc = p2.c - p1.c
    @assert 0.0 <= abs(dc) <= 1e6 "Distance (col coordinate = $(dc)) between $(p1) and $(p2) makes no sense "

    dh = p2.h - p1.h
    @assert 0.0 <= abs(dh) <= 1e6 "Distance (height coordinate = $(dh)) between $(p1) and $(p2) makes no sense "

    return sqrt(dr^2 + dc^2 + dh^2)
end


# """
# $(SIGNATURES)
# """
# function area(
#     p1::Tuple{Float64,Float64,Float64},
#     p2::Tuple{Float64,Float64,Float64},
#     p3::Tuple{Float64,Float64,Float64},
# )
#     a = dist(p1, p2)
#     b = dist(p2, p3)
#     c = dist(p3, p1)
#     s = (a + b + c) / 2.0

#     return sqrt(s * (s - a) * (s - b) * (s - c))
# end

"""
$(SIGNATURES)
"""
function area(v1::Vertex3D, v2::Vertex3D, v3::Vertex3D)
    a = dist(v1, v2)
    b = dist(v2, v3)
    c = dist(v3, v1)
    s = (a + b + c) / 2.0

    surface = s * (s - a) * (s - b) * (s - c)
    @assert surface >= 0.0 "Negative surface for $(v1),  $(v2),  $(v3)"

    return sqrt(s * (s - a) * (s - b) * (s - c))
end


"""
$(SIGNATURES)

A midpoint is the average between two points
"""
midpoint(p1::Vertex3D, p2::Vertex3D) =
    Vertex3D((p1.r + p2.r) / 2.0, (p1.c + p2.c) / 2.0, (p1.h + p2.h) / 2.0, 0.0, 0.0)


"""
$(SIGNATURES)

Centroid of three points.
"""
centroid(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D) = Vertex3D(
    (p1.r + p2.r + p3.r) / 3.0,
    (p1.c + p2.c + p3.c) / 3.0,
    (p1.h + p2.h + p3.h) / 3.0,
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

Velocity vector along `x` and `y` (velocity along h not necessary).

Implementation is a struc of arrays. Easier to vectorise.

"""
mutable struct FieldVertex3D
    size::Tuple{Int,Int}

    r::AbstractMatrix{Float64}
    c::AbstractMatrix{Float64}
    h::AbstractMatrix{Float64}

    vr::AbstractMatrix{Float64}
    vc::AbstractMatrix{Float64}

    rows_numbers::AbstractMatrix{Float64}
    cols_numbers::AbstractMatrix{Float64}

    # FieldVertex3D(size, mr, mc, mh, mvr, mvc, rows_numbers, cols_numbers) = new(size, mr, mc, mh, mvr, mvc, rows_numbers, cols_numbers)

end


function FieldVertex3D(height, width)
    # mr = zeros(Float64, height + 1, width + 1)
    # my = zeros(Float64, height + 1, width + 1)
    # mh = zeros(Float64, height + 1, width + 1)

    # mvr = zeros(Float64, height + 1, width + 1)
    # mvc = zeros(Float64, height + 1, width + 1)

    rows_numbers = repeat(Float64.(1:height+1), 1, width + 1)
    cols_numbers = repeat(Float64.(1:width+1)', height + 1, 1)

    mr = rows_numbers + rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
    mc = cols_numbers + rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
    mh = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
    mvr = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000
    mvc = rand(Float64, height + 1, width + 1) ./ 1_000 .- 0.5 / 1_000

    fv = FieldVertex3D((height, width), mr, mc, mh, mvr, mvc, rows_numbers, cols_numbers)
    reset_border_values!(fv)
    return fv
end


"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
function reset_border_values!(topleft::FieldVertex3D)
    # Reset the border at the fixed values fixed coordinates.
    topleft.r[1, :] .= topleft.rows_numbers[1, :]
    topleft.r[end, :] .= topleft.rows_numbers[end, :]
    topleft.r[:, 1] .= topleft.rows_numbers[:, 1]
    topleft.r[:, end] .= topleft.rows_numbers[:, end]

    topleft.c[1, :] .= topleft.cols_numbers[1, :]
    topleft.c[end, :] .= topleft.cols_numbers[end, :]
    topleft.c[:, 1] .= topleft.cols_numbers[:, 1]
    topleft.c[:, end] .= topleft.cols_numbers[:, end]

    fill_borders!(topleft.h, 0.0)
    fill_borders!(topleft.vr, 0.0)
    fill_borders!(topleft.vc, 0.0)
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
            fv.h[row, col],
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
function triangle3D(mesh::FaceMesh, row::Int, col::Int, side = Union{:top,:bottom})
    height, width = size(mesh)
    if 1 <= row <= height && 1 <= col <= width
        if side == :top
            t = mesh.toptriangles[row, col]
        elseif side == :bottom
            t = mesh.bottriangles[row, col]
        end

        t1_row, t1_col = t[1]
        t2_row, t2_col = t[2]
        t3_row, t3_col = t[3]

        r1 = mesh.topleft.r[t1_row, t1_col]
        c1 = mesh.topleft.c[t1_row, t1_col]
        h1 = mesh.topleft.h[t1_row, t1_col]
        vr1 = mesh.topleft.vr[t1_row, t1_col]
        vc1 = mesh.topleft.vc[t1_row, t1_col]
        @assert r1^2 + c1^2 + h1^2 <= 1e12 "Coordinate $(r1), $(c1), $(h1) of point #1 at $(row), $(col) side = $(side) makes no sense "

        r2 = mesh.topleft.r[t2_row, t2_col]
        c2 = mesh.topleft.c[t2_row, t2_col]
        h2 = mesh.topleft.h[t2_row, t2_col]
        vr2 = mesh.topleft.vr[t2_row, t2_col]
        vc2 = mesh.topleft.vc[t2_row, t2_col]
        @assert r2^2 + c2^2 + h2^2 <= 1e12 "Coordinate $(r2), $(c2), $(h2) of point #2 at $(row), $(col) side = $(side) makes no sense "

        r3 = mesh.topleft.r[t3_row, t3_col]
        c3 = mesh.topleft.c[t3_row, t3_col]
        h3 = mesh.topleft.h[t3_row, t3_col]
        vr3 = mesh.topleft.vr[t3_row, t3_col]
        vc3 = mesh.topleft.vc[t3_row, t3_col]
        @assert r3^2 + c3^2 + h3^2 <= 1e12 "Coordinate $(r3), $(c3), $(h3) of point #3 at $(row), $(col) side = $(side) makes no sense "

        p1 = Vertex3D(r1, c1, h1, vr1, vc1)
        p2 = Vertex3D(r2, c2, h2, vr2, vc2)
        p3 = Vertex3D(r3, c3, h3, vr3, vc3)
        return (p1, p2, p3)
    else
        # return (missing, missing, missing)
        @assert false "Coordinates out of bounds at $(row), $(col)"
    end
end

triangle3D(mesh::FaceMesh, ci::CartesianIndex{2}, side = Union{:top,:bottom}) =
    triangle3D(mesh, ci[1], ci[2], side)



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


max_sum_abs(m::AbstractMatrix) = (sum(abs.(m)) / length(m)) < 1_024
