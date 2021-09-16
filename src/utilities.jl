"""
$(TYPEDEF)

## Coordinates. Represent the data of each topleft corner of a pixel.
## posts/fences: One more corner than number of pixels.

"""
struct Vertex3D
    r::Float64
    c::Float64
    ϕ::Float64

    vr::Float64
    vc::Float64

    Vertex3D(r, c, ϕ, vr, vc) = new(r, c, ϕ, vr, vc)
    Vertex3D() = new(0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
$(SIGNATURES)
"""
function dist(p1::Tuple{Float64,Float64,Float64}, p2::Tuple{Float64,Float64,Float64})
    dr = p2[1] - p1[1]
    dc = p2[2] - p1[2]
    return sqrt(dr^2 + dc^2)
end

"""
$(SIGNATURES)
"""
function dist(p1::Vertex3D, p2::Vertex3D)
    dr = p2.r - p1.r
    dc = p2.c - p1.c
    return sqrt(dr^2 + dc^2)
end


"""
$(SIGNATURES)
"""
function area(v1::Vertex3D, v2::Vertex3D, v3::Vertex3D)
    a = dist(v1, v2)
    b = dist(v2, v3)
    c = dist(v3, v1)
    s = (a + b + c) / 2.0

    # surface_sq = (a + b + c) * (- a + b + c) * (a - b + c) * (a + b - c) / 16.0
    surface_sq = s * (s - a) * (s - b) * (s - c)
    return surface_sq <= 1e-100 ? 0.0 : sqrt(surface_sq)
end


"""
$(SIGNATURES)

Centroid of three points.
"""
centroid(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D) =
    Vertex3D((p1.r + p2.r + p3.r) / 3.0, (p1.c + p2.c + p3.c) / 3.0, 0.0, 0.0, 0.0)


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

Velocity vector along `row` and `col` (velocity along h not used).

Implementation is a struc of arrays. Easier to vectorise.

"""
mutable struct FieldVertex3D
    size::Tuple{Int,Int}

    r::AbstractMatrix{Float64}
    c::AbstractMatrix{Float64}

    # Velocity potential at each corner
    ϕ::AbstractMatrix{Float64}

    vr::AbstractMatrix{Float64}
    vc::AbstractMatrix{Float64}

    rows_numbers::AbstractMatrix{Float64}
    cols_numbers::AbstractMatrix{Float64}
end

# Information about each corner surrounding a pixel. Posts&fences warning!!!
function FieldVertex3D(height, width)
    rows_numbers = repeat(Float64.(1:height+1), 1, width + 1)
    cols_numbers = repeat(Float64.(1:width+1)', height + 1, 1)

    mr = rows_numbers
    mc = cols_numbers
    mϕ = zeros(Float64, height + 1, width + 1)

    mvr = zeros(Float64, height + 1, width + 1)
    mvc = zeros(Float64, height + 1, width + 1)

    fv = FieldVertex3D((height, width), mr, mc, mϕ, mvr, mvc, rows_numbers, cols_numbers)
    reset_border_values!(fv)
    return fv
end


"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
function reset_border_values!(corners::FieldVertex3D)
    # Reset the border at the fixed values fixed coordinates.
    corners.r[1, :] .= corners.rows_numbers[1, :]
    corners.r[end, :] .= corners.rows_numbers[end, :]

    corners.r[:, 1] .= corners.rows_numbers[:, 1]
    corners.r[:, end] .= corners.rows_numbers[:, end]

    corners.c[1, :] .= corners.cols_numbers[1, :]
    corners.c[end, :] .= corners.cols_numbers[end, :]
    corners.c[:, 1] .= corners.cols_numbers[:, 1]
    corners.c[:, end] .= corners.cols_numbers[:, end]

    fill_borders!(corners.vr, 0.0)
    fill_borders!(corners.vc, 0.0)
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
            fv.ϕ[row, col],
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
    corners::FieldVertex3D

    toptriangles::AbstractMatrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}
    bottriangles::AbstractMatrix{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}}}

    FaceMesh(s, tl, tt, bt) = new(s, tl, tt, bt)

    function FaceMesh(height::Int, width::Int)
        m_corner = FieldVertex3D(height, width)
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

        return FaceMesh((height, width), m_corner, t_triangles, b_triangles)
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
function triangle3D(mesh::FaceMesh, row::Int, col::Int, side = Union{:top,:bottom})
    height, width = size(mesh)
    max_distance_squared = (height + 1)^2 + (width + 1)^2

    if 1 <= row <= height && 1 <= col <= width
        if side == :top
            t = mesh.toptriangles[row, col]
        elseif side == :bottom
            t = mesh.bottriangles[row, col]
        end

        t1_row, t1_col = t[1]
        t2_row, t2_col = t[2]
        t3_row, t3_col = t[3]

        r1 = mesh.corners.r[t1_row, t1_col]
        c1 = mesh.corners.c[t1_row, t1_col]
        ϕ1 = mesh.corners.ϕ[t1_row, t1_col]
        vr1 = mesh.corners.vr[t1_row, t1_col]
        vc1 = mesh.corners.vc[t1_row, t1_col]
        @assert r1^2 + c1^2 <= max_distance_squared "Coordinate $(r1), $(c1), $(ϕ1) of point #1 at $(row), $(col), side = $(side) makes no sense "

        r2 = mesh.corners.r[t2_row, t2_col]
        c2 = mesh.corners.c[t2_row, t2_col]
        ϕ2 = mesh.corners.ϕ[t2_row, t2_col]
        vr2 = mesh.corners.vr[t2_row, t2_col]
        vc2 = mesh.corners.vc[t2_row, t2_col]
        @assert r2^2 + c2^2 <= max_distance_squared "Coordinate $(r2), $(c2), $(ϕ2) of point #2 at $(row), $(col), side = $(side) makes no sense "

        r3 = mesh.corners.r[t3_row, t3_col]
        c3 = mesh.corners.c[t3_row, t3_col]
        ϕ3 = mesh.corners.ϕ[t3_row, t3_col]
        vr3 = mesh.corners.vr[t3_row, t3_col]
        vc3 = mesh.corners.vc[t3_row, t3_col]
        @assert r3^2 + c3^2 <= max_distance_squared "Coordinate $(r3), $(c3), $(ϕ3) of point #3 at $(row), $(col), side = $(side) makes no sense "

        p1 = Vertex3D(r1, c1, ϕ1, vr1, vc1)
        p2 = Vertex3D(r2, c2, ϕ2, vr2, vc2)
        p3 = Vertex3D(r3, c3, ϕ3, vr3, vc3)
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

"""
$(SIGNATURES)
"""
area(t::Tuple{Vertex3D,Vertex3D,Vertex3D}) = area(t...)

"""
$(SIGNATURES)

A Mesh is a collection of triangles. The brightness flowing through a given triangle is just proportional to its
area in the r, r plane. h is ignored.

The function returns a matrix with the quantity of light coming from each 'rectangle'  around a corner. That 'rectangle'
has been shifted and flexed around.
"""
function get_lens_pixels_area(mesh::FaceMesh)
    height, width = size(mesh)

    top_tri_area =
        area.([triangle3D(mesh, row, col, :bottom) for row ∈ 1:height, col ∈ 1:width])
    bot_tri_area =
        area.([triangle3D(mesh, row, col, :bottom) for row ∈ 1:height, col ∈ 1:width])

    total_area = top_tri_area + bot_tri_area

    average_energy_per_pixel = average(total_area)
    return total_area / average_energy_per_pixel
end
