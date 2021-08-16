using Images


# This implements the method of caustics control described in this paper:
# https://www.researchgate.net/profile/Yonghao_Yue/publication/274483217_Poisson-Based_Continuous_Surface_Generation_for_Goal-Based_Caustics/links/575b4ceb08ae414b8e467a5f.pdf


mutable struct Point3D
    x::Float64
    y::Float64
    z::Float64
    ix::Int
    iy::Int
end

struct Triangle
    pt1::Int64
    pt2::Int64
    pt3::Int64
end

struct Mesh
    nodes::Vector{Point3D}
    nodeArray::Matrix{Point3D}
    triangles::Vector{Triangle}
    width::Int
    height::Int
end

function squareMesh(width::Int, height::Int)
    # This func returns a square mesh, centered on zero, with (width * height) nodes
    nodeList = Vector{Point3D}(undef, height * width)
    nodeArray = Matrix{Point3D}(undef, width, height)
    count = 1
    midpoint = width / 2
    for y = 1:height
        for x = 1:width
            newPoint = Point3D(x, y, 0, x, y)
            nodeList[count] = newPoint
            nodeArray[x, y] = newPoint
            count += 1
        end
    end

    triangles = Vector{Triangle}(undef, (width - 1) * (height - 1) * 2)
    count = 1
    for y = 1:(height - 1)
        for x = 1:(width - 1)
          # here x and y establish the column of squares we're in
            index_ul = (y - 1) * width + x
            index_ur = index_ul + 1

            index_ll = y * width + x
            index_lr = index_ll + 1

            triangles[count] = Triangle(index_ul, index_ll, index_ur)
            count += 1
            triangles[count] = Triangle(index_lr, index_ur, index_ll)
            count += 1
        end
    end

    newMesh = Mesh(nodeList, nodeArray, triangles, width, height)
end

function dist(p1::Point3D, p2::Point3D)
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    sqrt(dx * dx + dy * dy)
end

function midpoint(p1::Point3D, p2::Point3D)
  # a midpoint is just the average between two points
    Point3D(.5p1.x + .5p2.x, .5p1.y + .5p2.y, .5p1.z + .5p2.z, 0, 0)
end

function centroid(mesh::Mesh, index::Int)
  # Warning: not guaranteed to work in 3D?
    triangle = mesh.triangles[index]
    p1 = mesh.nodes[triangle.pt1]
    p2 = mesh.nodes[triangle.pt2]
    p3 = mesh.nodes[triangle.pt3]
    centroid(p1, p2, p3)
end

function centroid(p1::Point3D, p2::Point3D, p3::Point3D)
    Point3D(1 / 3 * (p1.x + p2.x + p3.x), 1 / 3 * (p1.y + p2.y + p3.y), 1 / 3 * (p1.z + p2.z + p3.z), 0, 0)
end

function findT(p1::Point3D, p2::Point3D, p3::Point3D, dp1::Point3D, dp2::Point3D, dp3::Point3D)
    # Given 3 points and 3 velocities, calculate the t required to bring the area of that triangle to zero
    x1 = p2.x - p1.x
    y1 = p2.y - p1.y

    x2 = p3.x - p1.x
    y2 = p3.y - p1.y

    u1 = dp2.x - dp1.x
    v1 = dp2.y - dp1.y

    u2 = dp3.x - dp1.x
    v2 = dp3.y - dp1.y

    a = u1 * v2 - u2 * v1
    b = x1 * v1 + y2 * u1 - x2 * v1 - y1 * u2
    c = x1 * y2 - x2 * y1
    if a != 0
        quotient = b^2 - 4a * c
        if quotient >= 0
            d = sqrt(quotient)
            (-b - d) / 2a, (-b + d) / 2a
        else
            -123.0, -123.0
        end
    else
        # cool, there just isn't any dependence on t^2, but there is still on t!
        -c / b, -c / b
    end
end

function triangle_area(mesh::Mesh, index::Int)
    triangle = mesh.triangles[index]
    pt1 = mesh.nodes[triangle.pt1]
    pt2 = mesh.nodes[triangle.pt2]
    pt3 = mesh.nodes[triangle.pt3]
    triangle_area(pt1, pt2, pt3)
end

function triangle_area(p1::Point3D, p2::Point3D, p3::Point3D)
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a + b + c) / 2
    sqrt(s * (s - a) * (s - b) * (s - c))
end


function saveObj(mesh::Mesh, filename::String, scale=1.0, scalez=1.0)
  # This function saves the mesh object in stl format
    open(filename, "w") do io
        for vertex in mesh.nodes
            println(io, "v ", vertex.x * scale, " ", vertex.y * scale, " ", vertex.z * scalez)
        end

        for face in mesh.triangles
            println(io, "f ", face.pt1, " ", face.pt2, " ", face.pt3)
        end

        println(io, "dims ", mesh.width, " ", mesh.height)
    end
end

function Obj2Mesh(filename)
    lines = readlines(filename)

    vertexLines = [l for l in lines if startswith(l, "v")]
    nodeList = Vector{Point3D}(undef, size(vertexLines))
    count = 1
    for line in vertexLines
        elements = split(line, " ")
        x = parse(Float64, elements[2])
        y = parse(Float64, elements[3])
        z = parse(Float64, elements[4]) * 10
        pt = Point3D(x, y, z, 0, 0)
        nodeList[count] = pt
        count += 1
    end

    faceLines = [l for l in lines if startswith(l, "f")]
    triangles = Vector{Triangle}(undef, size(faceLines))
    for line in faceLines
        elements = split(line, " ")
        triangle = Triangle(parse(Int64, elements[2]), parse(Int64, elements[3]), parse(Int64, elements[4]))
    end

    dimsLines = [l for l in lines if startswith(l, "dims")]
    elements = split(dimsLines[1], " ")

    newMesh = Mesh(nodeList, triangles, parse(Int64, elements[2]), parse(Int64, elements[3]))
end

function ∇(f::Matrix{Float64})
    w, h = size(f)
    ∇fᵤ = Matrix{Float64}(undef, w, h)   # the right edge will be filled with zeros
    ∇fᵥ = Matrix{Float64}(undef, w, h)   # the buttom edge will be filled with zeros

    for x = 1:w
        for y = 1:h
            if x == w
                ∇fᵤ[x, y] = 0
            else
                ∇fᵤ[x, y] = f[x + 1, y] - f[x, y]
            end
        end
    end

    for x = 1:w
        for y = 1:h
            if y == h
                ∇fᵥ[x, y] = 0
            else
                ∇fᵥ[x, y] = f[x, y + 1] - f[x, y]
            end
        end
    end

    return ∇fᵤ, ∇fᵥ
end



function getPixelArea(mesh::Mesh)
    # A Mesh is a grid of 3D points. The X and Y coordinates are not necessarily aligned or square
    # The Z coordinate represents the value. brightness is just proportional to area.
    pixelAreas = Matrix{Float64}(undef, mesh.width-1, mesh.height-1)
    for x = 1:mesh.width-1
        for y = 1:mesh.height-1
            upperLeft = mesh.nodeArray[x, y]
            upperRight = mesh.nodeArray[x + 1, y]
            lowerLeft = mesh.nodeArray[x, y + 1]
            lowerRight = mesh.nodeArray[x + 1, y + 1]

            #=
            *------*
            |    / |
            |   /  |
            |  /   |
            | /    |
            *------*
            =#
            area = triangle_area(lowerLeft, upperRight, upperLeft) + triangle_area(lowerLeft, lowerRight, upperRight)

            pixelAreas[x, y] = area
        end
    end
    pixelAreas
end


function relax!(matrix::Matrix{Float64}, D::Matrix{Float64})
    # This function implements successive over relaxation for a matrix and its associated error matrix
    # There is a hardcoded assumption of Neumann boundary conditions--that the derivative across the
    # boundary must be zero in all cases. See:
    # https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions
    # sz = size(matrix)
    # width = sz[1]
    # height = sz[2]
    width, height = size(matrix)
    # ω = 2 / (1 + π / width)
    ω = 1.99
    # println("OMEGA $(ω)")
    max_update = 0
    for y = 1:height
        for x = 1:width
            val = matrix[x, y]

            if x == 1 && y == 1
                # Top left corner
                val_down = matrix[x, y + 1]
                val_right = matrix[x + 1, y]
                delta = ω / 2 * (val_down + val_right - 2 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif x == 1 && y == height
                # Bottom left corner
                val_up = matrix[x, y - 1]
                val_right = matrix[x + 1, y]
                delta = ω / 2 * (val_up + val_right - 2 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif x == width && y == 1
                # Top right corner
                val_down = matrix[x, y + 1]
                val_left = matrix[x - 1, y]
                delta = ω / 2 * (val_down + val_left - 2 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif x == width && y == height
                # Bottom right corner
                val_up = matrix[x, y - 1]
                val_left = matrix[x - 1, y]
                delta = ω / 2 * (val_up + val_left - 2 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta

            elseif x == 1
                # Along the left edge, but not the top or buttom corner
                val_up = matrix[x, y - 1]
                val_down = matrix[x, y + 1]
                val_right = matrix[x + 1, y]
                delta = ω / 3 * (val_up + val_down + val_right - 3 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif x == width
                # Along the right edge, but not the top or buttom corner
                val_up = matrix[x, y - 1]
                val_down = matrix[x, y + 1]
                val_left = matrix[x - 1, y]
                delta = ω / 3 * (val_up + val_down + val_left - 3 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif y == 1
                # Along the top edge, but not the left or right corner
                val_down = matrix[x, y + 1]
                val_left = matrix[x - 1, y]
                val_right = matrix[x + 1, y]
                delta = ω / 3 * (val_down + val_left + val_right - 3 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            elseif y == height
                # Along the bottom edge, but not the left or right corner
                val_up = matrix[x, y - 1]
                val_left = matrix[x - 1, y]
                val_right = matrix[x + 1, y]
                delta = ω / 3 * (val_up + val_left + val_right - 3 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            else
                # The normal case, in the middle of the mesh!
                val_up = matrix[x, y - 1]
                val_down = matrix[x, y + 1]
                val_left = matrix[x - 1, y]
                val_right = matrix[x + 1, y]

                # The new way
                # ∇x₁ =


                # The old way
                delta = ω / 4 * (val_up + val_down + val_left + val_right - 4 * val - D[x, y])
                if abs(delta) > max_update
                    max_update = abs(delta)
                end
                matrix[x, y] += delta
            end
            # node.z = .25 * (node_up.z + node_down.z + node_left.z + node_right.z) # simple averaging
            # node.z += ω/4 * (node_up.z + node_down.z + node_left.z + node_right.z - 4 * node.z)

            # matrix[x, y] += ω/4 * (val_up + val_down + val_left + val_right - 4 * val - D[x, y])
        end
    end

    max_update

    # for y = 1:height
    #     for x = 1:width
    #         val = matrix[x, y]
    #     end
    # end
end

function matrix_to_mesh(matrix::Matrix{Float64})
    # This function takes a 512x512 matrix and returns a 512x512 mesh
    w, h = size(matrix)
    retval = squareMesh(w, h)
    for x = 1:w
        for y = 1:h
            index = (y - 1) * retval.width + x
            node = retval.nodes[index]
            node.z = matrix[x, y]
        end
    end
    retval
end



function marchMesh!(mesh::Mesh, ϕ::Matrix{Float64})
    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)

    imgWidth, imgHeight = size(ϕ)   # should be 512x512

    # For each point in the mesh we need to figure out its velocity
    velocities = Matrix{Point3D}(undef, mesh.width, mesh.height)

    for x in 1:mesh.width
        for y in 1:mesh.height
            # XY coordinates in the mesh ARE XY coordinates in the image. The mesh just needs an extra row and column
            # at the bottom right edge so that the triangles can be closed


            if x == mesh.width
                u = 0
            else
                if y == mesh.height
                    u = ∇ϕᵤ[x, y - 1]
                else
                    u = ∇ϕᵤ[x, y]
                end
            end

            if y == mesh.height
                v = 0
            else
                if x == mesh.width
                    v = ∇ϕᵥ[x - 1, y]
                else
                    v = ∇ϕᵥ[x, y]
                end
            end
            velocities[x, y] = Point3D(-u, -v, 0, 0, 0)
        end
    end


    min_t = 10000
    triangleCount = 1
    for triangle in mesh.triangles
        p1 = mesh.nodes[triangle.pt1]
        p2 = mesh.nodes[triangle.pt2]
        p3 = mesh.nodes[triangle.pt3]

        v1 = velocities[p1.ix, p1.iy]
        v2 = velocities[p2.ix, p2.iy]
        v3 = velocities[p3.ix, p3.iy]

        t1, t2 = findT(p1, p2, p3, v1, v2, v3)

        if t1 > 0 && t1 < min_t
            min_t = t1
        end

        if t2 > 0 && t2 < min_t
            min_t = t2
        end
        triangleCount += 1
    end

    println("Overall min_t:", min_t)
    δ = min_t / 2

    for point in mesh.nodes
        v = velocities[point.ix, point.iy]
        point.x = v.x * δ + point.x
        point.y = v.y * δ + point.y
    end

    # saveObj(mesh, "gateau.obj")
end

function quantifyLoss(D, suffix, img)
    println("Loss:")
    println("Minimum: $(minimum(D))")
    println("Maximum: $(maximum(D))")
    blue = zeros(size(D))
    blue[D .> 0] = D[D .> 0]

    red = zeros(size(D))
    red[D .< 0] = -D[D .< 0]
    green = zeros(size(D))

    println(size(blue))
    println(size(red))
    println(size(green))
    rgbImg = RGB.(red, green, blue)'
    save("loss_$(suffix).png", map(clamp01nan, rgbImg))

    # println("Saving output image:")
    # println(typeof(img))
    # E = Gray.(D)
    # println(typeof(E))
    # outputImg = img - E
    # save("actual_$(suffix).png", outputImg)
end

function oneIteration(meshy, img, suffix)
    # remember meshy is 512x512 just like the image 512x512
    # so LJ is 512x512
    LJ = getPixelArea(meshy)
    D = Float64.(LJ - img)
    
    # Save the loss image as a png
    println(minimum(D))
    println(maximum(D))
    quantifyLoss(D, suffix, img)
    # save("loss_$(suffix).png", colorview(Gray, D))
    # return
    width, height = size(img)

    ϕ = Matrix{Float64}(undef, width, height)

    for i = 1:10240/2
        max_update = relax!(ϕ, D)

        if i % 500 == 0
            println(max_update)
        end
        if max_update < 0.00001
            println("Convergence reached at step $(i) with max_update of $(max_update)")
            break
        end
    end

    saveObj(matrix_to_mesh(ϕ * .02), "phi_$(suffix).obj")
    # saveObj(matrix_to_mesh(D * 10), "D_$(suffix).obj")

    # Now we need to march the x,y locations in our mesh according to this gradient!
    marchMesh!(meshy, ϕ)
end

function setHeights!(mesh, heights, heightScale=1.0, heightOffset=50)
    width, height = size(heights)
    for y = 1:height
        for x = 1:width
            mesh.nodeArray[x, y].z = heights[x, y] * heightScale + heightOffset
            if x == 100 && y == 100
                println("Example heights: $(heights[x, y])  and  $(heights[x, y] * heightScale) and $(heights[x, y] * heightScale + heightOffset)")
            end
        end
    end

    # get the side edge
    for y = 1:height
        mesh.nodeArray[width+1, y].z = mesh.nodeArray[width, y].z
    end

    # get the bottom edge
    for x = 1:width+1
        mesh.nodeArray[x, height+1].z = mesh.nodeArray[x, height].z
    end


    # # get the pesky corner!
    # mesh.nodeArray[width + 1, height + 1].z = mesh.nodeArray[width, height].z
end


function setHeights(mesh, heights)
    width = mesh.width
    height = mesh.height
    nodes = Vector{Point3D}(undef, size(mesh.nodes))
    nodeArray = Matrix{Point3D}(undef, 0, 0)
    triangles = Vector{Triangle}(undef, size(mesh.triangles))
    scale = 1
    count = 1
    w, h = size(heights)
    for y = 1: height
        for x = 1: width
            count += 1
            point = mesh.nodes[count]

            z = heights[x, y]

            newPoint = Point3D(point.x * scale, point.y * scale, z * scale, 0, 0)
            nodes[count] = newPoint
        end
    end

    Mesh(nodes, nodeArray, triangles, width, height)
end

function solidify(inputMesh, offset=100)
    width = inputMesh.width
    height = inputMesh.height
    totalNodes = width * height * 2
    nodeList = Vector{Point3D}(undef, totalNodes)
    nodeArrayTop = Matrix{Point3D}(undef, width, height)
    nodeArrayBottom = Matrix{Point3D}(undef, width, height)

    # imagine a 4x4 image. 4 * 2 + 2 * 2 = 12
    numEdgeNodes = width * 2 + (height - 2) * 2

    numTrianglesTop = (width-1)*(height-1) * 2
    numTrianglesBottom = numTrianglesTop
    numTrianglesEdges = numEdgeNodes * 2

    totalTriangles = numTrianglesBottom + numTrianglesTop + numTrianglesEdges

    println("Specs: $(width)  $(height)  $(totalNodes)  $(numEdgeNodes)  $(numTrianglesBottom) $(totalTriangles)")

    # Build the bottom surface
    count = 1
    for y = 1:height
        for x = 1:width
            newPoint = Point3D(x, y, -offset, x, y)
            nodeList[count] = newPoint
            nodeArrayBottom[x, y] = newPoint
            count += 1
        end
    end

    # Copy in the top surface
    for y = 1:height
        for x = 1:width
            node = inputMesh.nodeArray[x, y]
            copiedPoint = Point3D(node.x, node.y, node.z, node.ix, node.iy)
            if node.ix != x
                println("OH NO POINTS NOT MATCHED $(x) vs $(node.ix)")
            end
            if node.iy != y
                println("OH NO POINTS NOT MATCHED $(y) vs $(node.iy)")
            end

            nodeList[count] = copiedPoint
            nodeArrayTop[x, y] = copiedPoint
            count += 1
        end
    end

    println("We now have $(count-1) valid nodes")

    triangles = Vector{Triangle}(undef, totalTriangles)
    # Build the triangles for the bottom surface
    count = 1
    for y = 1:(height - 1)
        for x = 1:(width - 1)
          # here x and y establish the column of squares we're in
            index_ul = (y - 1) * width + x
            index_ur = index_ul + 1

            index_ll = y * width + x
            index_lr = index_ll + 1

            triangles[count] = Triangle(index_ul, index_ll, index_ur)
            count += 1
            triangles[count] = Triangle(index_lr, index_ur, index_ll)
            count += 1
        end
    end

    println("We've filled up $(count-1) triangles")
    if count != numTrianglesBottom + 1
        println("Hmm aren't count and triangles bottom equal? $(count) vs $(numTrianglesBottom + 1)")
    end

    # Build the triangles for the top surface
    for y = 1:(height - 1)
        for x = 1:(width - 1)
          # here x and y establish the column of squares we're in
            index_ul = (y - 1) * width + x + totalNodes / 2
            index_ur = index_ul + 1

            index_ll = y * width + x + totalNodes / 2
            index_lr = index_ll + 1

            triangles[count] = Triangle(index_ul, index_ur, index_ll)
            count += 1
            triangles[count] = Triangle(index_lr, index_ll, index_ur)
            count += 1
        end
    end

    println("We've filled up $(count-1) triangles")

    # Build the triangles to close the mesh
    x = 1
    for y = 1:(height - 1)
        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = y * width + x
        ur = lr + totalNodes / 2
        triangles[count] = Triangle(ll, ul, ur)
        count += 1
        triangles[count] = Triangle(ur, lr, ll)
        count += 1
    end

    x = width
    for y = 1:(height - 1)
        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = y * width + x
        ur = lr + totalNodes / 2
        triangles[count] = Triangle(ll, ur, ul)
        count += 1
        triangles[count] = Triangle(ur, ll, lr)
        count += 1
    end

    y = 1
    for x = 2: width
        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = (y - 1) * width + (x - 1)
        ur = lr + totalNodes / 2
        triangles[count] = Triangle(ll, ul, ur)
        count += 1
        triangles[count] = Triangle(ur, lr, ll)
        count += 1
    end

    y = height
    for x = 2: width
        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = (y - 1) * width + (x - 1)
        ur = lr + totalNodes / 2
        triangles[count] = Triangle(ll, ur, ul)
        count += 1
        triangles[count] = Triangle(ur, ll, lr)
        count += 1
    end

    Mesh(nodeList, nodeArrayBottom, triangles, width, height)
end

function findSurface(mesh, image, f, imgWidth)
    width, height = size(image)

    # imgWidth = .1 # m
    # f = 1.0  # m
    H = f
    metersPerPixel = imgWidth / width
    println(metersPerPixel)

    # η = 1.49
    n₂ = 1
    n₁ = 1.49
    Nx = Matrix{Float64}(undef, width + 1, height + 1)
    Ny = Matrix{Float64}(undef, width + 1, height + 1)

    for j = 1:height
        for i = 1:width
            node = mesh.nodeArray[i, j]
            dx = (node.ix - node.x) * metersPerPixel
            dy = (node.iy - node.y) * metersPerPixel

            little_h = node.z * metersPerPixel
            H_minus_h = H - little_h
            dz = H_minus_h


            # k = η * sqrt(dx * dx + dy * dy + H_minus_h * H_minus_h) - H_minus_h
            # Nx[i, j] = 1/k * dx
            # Ny[i, j] = 1/k * dy
            Ny[i, j] = tan(atan(dy / dz) / (n₁ - 1))
            Nx[i, j] = tan(atan(dx / dz) / (n₁ - 1))


        end
    end

    divergence = Matrix{Float64}(undef, width, height)
    # We need to find the divergence of the Vector field described by Nx and Ny

    for j = 1:height
        for i = 1:width
            δx = (Nx[i+1, j] - Nx[i, j])
            δy = (Ny[i, j+1] - Ny[i, j])
            divergence[i, j] = δx + δy

            if i == 100 && j == 100
                println("div: $(divergence[i, j])")
            end
        end
    end
    println("Have all the divergences")

    h = Matrix{Float64}(undef, width, height)
    max_update = 0
    for i = 1:10240/4
        max_update = relax!(h, divergence)

        if i % 100 == 0
            println(max_update)
        end
        if max_update < 0.00001
            println("Convergence reached at step $(i) with max_update of $(max_update)")
            break
        end
    end
    # saveObj(matrix_to_mesh(h / 10), "heightmap.obj")
    h, metersPerPixel
end

function testSquareMesh()
    mesh = squareMesh(100, 50)

    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodeArray[1, 1].x = 8
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodes[1].y += 12
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

end

function testSolidify()
    println("Testing solidification")
    width = 100
    height = 100
    origMesh = squareMesh(width, height)

    for y = 1: height
        for x = 1: width
            x2 = (x - width/2) / width
            y2 = (y - height/2) / height
            value = x2 * x2 + y2 * y2
            origMesh.nodeArray[x, y].z = 15 - value * 25
        end
    end

    saveObj(origMesh, "testSolidify.obj")
    solidMesh = solidify(origMesh, 0)
    saveObj(solidMesh, "testSolidify2.obj")
end

function main()
    # img = Gray.(load("cat.jpg"))
    # img = Gray.(load("necco2.jpg"))
    img = Gray.(load("cat_posing.jpg"))
    img2 = permutedims(img) * 1.0
    width, height = size(img2)

    # meshy is the same size as the image
    meshy = squareMesh(width + 1, height + 1)
    # We need to boost the brightness of the image so that its sum and the sum of the area are equal
    mesh_sum = width * height
    image_sum = sum(img2)
    boost_ratio = mesh_sum / image_sum
    # img3 is 512x512
    img3 = img2 .* boost_ratio    
    oneIteration(meshy, img3, "it1")
    oneIteration(meshy, img3, "it2")
    oneIteration(meshy, img3, "it3")
    # oneIteration(meshy, img3, "it4")
    # oneIteration(meshy, img3, "it5")
    # oneIteration(meshy, img3, "it6")

    artifactSize = 0.15  # meters

    h, metersPerPixel = findSurface(meshy, img3, 1.0, artifactSize)

    setHeights!(meshy, h, 1.0)
    # newMesh = setHeights(meshy, h)

    solidMesh = solidify(meshy)
    saveObj(solidMesh, "final_scripted.obj", 1/512.0 * artifactSize, 1/512.0 * artifactSize)
    meshy, img3
end


main()

