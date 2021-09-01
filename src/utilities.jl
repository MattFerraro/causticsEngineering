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

    Point3D() = new(0.0, 0.0, 0.0, 0.0, 0.0)
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
$(TYPEDEF)


Meshes representing the block of acrylate. Bottom is facing the light source; Top is facing the caustics.
CHECK whether bottomMesh is at all useful.
"""
struct TopBottomMeshes
    height::Int
    width::Int

    bottomMesh::Mesh
    topMesh::Mesh

    TopBottomMeshes(height::Int, width::Int) =
        new(height, width, Mesh(height, width), Mesh(height, width))
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

"""
$(SIGNATURES)
"""
function centroid(mesh::Mesh, height::Int, width::Int)
    mesh_height, _ = size(mesh.topTriangles)
    index = height + width * mesh_height

    return centroid(mesh, index)
end
