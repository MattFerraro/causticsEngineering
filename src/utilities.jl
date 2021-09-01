"""
$(TYPEDEF)
"""
mutable struct Point3D
    # Coordinates
    x::Float64
    y::Float64
    z::Float64

    # Velocity vector (velocity along z not necessary)
    vx::Float64
    vy::Float64

    Point3D() = new(0.0, 0.0, Height_Offset, 0.0, 0.0)
end

"""
$(SIGNATURES)
"""
function dist(p1::Point3D, p2::Point3D)
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dz = p2.z - p1.z
    return sqrt(dx^2 + dy^2 + dz^2)
end


"""
$(SIGNATURES)

A midpoint is the average between two points
"""
midpoint(p1::Point3D, p2::Point3D) =
    Point3D((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0, (p1.z + p2.z) / 2.0, 0.0, 0.0)


"""
$(SIGNATURES)

Barycentre of three points.
"""
centroid(p1::Point3D, p2::Point3D, p3::Point3D) = Point3D(
    (p1.x + p2.x + p3.x) / 3.0,
    (p1.y + p2.y + p3.y) / 3.0,
    (p1.z + p2.z + p3.z) / 3.0,
    0,
    0,
)


"""
$(SIGNATURES)
"""
function triangle_area(p1::Point3D, p2::Point3D, p3::Point3D)
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

`
*------*
|    / |
|   /  |
|  /   |
| /    |
*------*
`


Although the corners of each triangle match corners of the rectangles, it is sometimes easier to think in terms
of triangles rather than corners coming from different rectangles.
"""
struct Triangle
    pt1::Point3D
    pt2::Point3D
    pt3::Point3D

    Triangle() = new(Point3D(), Point3D(), Point3D())
    Triangle(pt1, pt2, pt3) = new(pt1, pt2, pt3)
end

"""
$(SIGNATURES)
"""
centroid(t::Triangle) = centroid(t.p1, t.p2, t.p3)

"""
$(SIGNATURES)
"""
triangle_area(t::Triangle) = triangle_area(t.p1, t.p2, t.p3)



"""
$(TYPEDEF)

A mesh is represents a single surface. First the surface is split into rectangles of size (width, height). Then,
each rectangle is broken into a top-left triangle and a bottom-right triangle.

`rectangles` represents the matrix of rectangles. Since each corner is shared by 4 rectangles, it is enough to
only record the top-left corner of each.
"""
struct Mesh
    rectangles::Matrix{Point3D}

    Mesh(height::Int, width::Int) = new(Matrix{Point3D}(undef, height + 1, width + 1))
end



"""
$(SIGNATURES)

Return the dimension of the mesh as the number of rectangles. It does not return the number of corner points.
"""
function base.size(mesh::Mesh)
    height, width = size(mesh.rectangles)
    return height - 1, width - 1
end

"""
$(SIGNATURES)

Warning: not guaranteed to work in 3D?
"""
function centroid(mesh::Mesh, index::Int)
    triangle = mesh.topTriangles[index]
    p1 = mesh.topNodes[triangle.pt1]
    p2 = mesh.topNodes[triangle.pt2]
    p3 = mesh.topNodes[triangle.pt3]

    return centroid(p1, p2, p3)
end


function top_triangle(mesh::Mesh, row::Int, col::Int)
    height, width = size(mesh)

    if (1 <= row <= height) && (1 <= col <= width)
        # The order is important
        pt1 = mesh.rectangles[row, col]
        pt2 = mesh.rectangles[row+1, col]
        pt3 = mesh.rectangles[row, col+1]
        return Triangle(pt1, pt2, pt3)
    else
        return nothing
    end
end

function bottom_triangle(mesh::Mesh, row::Int, col::Int)
    height, width = size(mesh)

    if (1 <= row <= height) && (1 <= col <= width)
        # The order is important
        pt1 = mesh.rectangles[row, col+1]
        pt2 = mesh.rectangles[row+1, col]
        pt3 = mesh.rectangles[row+1, col+1]
        return Triangle(pt1, pt2, pt3)
    else
        return nothing
    end
end
