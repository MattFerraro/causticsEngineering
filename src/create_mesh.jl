
"""
$(SIGNATURES)

This func returns a square mesh, centered on zero, with (width * height) nodes
"""
function create_mesh(width::Int, height::Int)

    new_mesh = Mesh(height, width)

    for row = 1:height, col = 1:width
        new_mesh.rectangles[row, col] = Point3D(Float64(row), Float64(col), 0.0, 0.0, 0.0)

        new_mesh.topTriangles[row, col] = Triangle(
            Point3D(Float64(row), Float64(col), 0.0, 0.0, 0.0),
            Point3D(Float64(row+1), Float64(col), 0.0, 0.0, 0.0),
            Point3D(Float64(row), Float64(col+1), 0.0, 0.0, 0.0))
        new_mesh.topNodes[row, col] = centroid(new_mesh.topTriangles[row, col])

        new_mesh.botTriangles[row, col] = Triangle(
            Point3D(Float64(row), Float64(col+1), 0.0, 0.0, 0.0),
            Point3D(Float64(row+1), Float64(col), 0.0, 0.0, 0.0),
            Point3D(Float64(row+1), Float64(col+1), 0.0, 0.0, 0.0))
        new_mesh.botNodes[row, col] = centroid(new_mesh.botTriangles[row, col])
    end

    return new_mesh
end


"""
$(SIGNATURES)

Given 3 points and 3 velocities, calculate the `t` required to bring the area of that triangle to zero
"""
function find_maximum_t(
    p1::Point3D,
    p2::Point3D,
    p3::Point3D,
)

    # The area of a triangle is 1/2 * [ Ax (By - Cy) + Bx (Cy - Ay) + Cx (Ay - By)]
    # where each point of the triangle is where it will be after time t
    # i.e. a point goes from P to P+tV where V is the velocity of that point.
    # Notation is here B -> B + t_vB

    # To make the calculation simpler, everything is translated so that A is at the
    # origin of the plane and its velocity is nil.
    Bx = p2.x - p1.x; By = p2.y - p1.y
    Cx = p3.x - p1.x; Cy = p3.y - p1.y

    t_vBx = p2.vx - p1.vx; t_vBy = p2.vy - p1.vy
    t_vCx = p3.vx - p1.vx; t_vCy = p3.vy - p1.vy

    # After this, given that Ax = Ay = 0, the area is nil iff
    # Bx Cy - Cx By = 0.
    # After expansion and reshuffling to have a quadratic equation where t
    # is the variable, the coefficients of that equation are:

    a = t_vBy * t_vCx - t_vCy * t_vBx
    b = By * t_vCx + Cx * t_vBy - (Bx * t_vCy + Cy * t_vBx)
    c = Cx * By - Bx * Cy

    # if a = 0, this is just a linear equation.
    if a == 0
        # Fingers crossed that b != 0
        return -c / b, -c / b
    else

        discriminant = b^2 - 4a * c

        # If there is a solution
        if discriminant >= 0
            d = sqrt(discriminant)
            return (-b - d) / 2a, (-b + d) / 2a
        else
            # There can be no solution if, after translation, B abd C move in parallel direction.
            # C will never end up on the line AB.result_mesh.rectangles[row, col].z = height_map[row, col]
        result_mesh.rectangles[row, col].z = height_map[row, col]

            # Very unlikely with Float64.
            return missing, missing
        end
    end
end


"""
$(SIGNATURES)

This function saves the mesh object in stl format.
"""
function save_stl!(
    mesh::Mesh,
    filename::String;
    scale = 1.0,
    scalez = 1.0,
    reverse = false,
    flipxy = false,
)

    height, width = size(mesh.rectangles)

    open(filename, "w") do io
        for row ∈ 1:height, col ∈ 1:width
            n = mesh.topNodes[row, col]
            if flipxy
                println(io,"v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io,"v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            t = mesh.topTriangles[row, col]
            if reverse
                println(io, "f $(t.pt3) $(t.pt2) $(t.pt1)")
            else
                println(io, "f $(t.pt1) $(t.pt2) $(t.pt3)")
            end

            n = mesh.botNodes[row, col]
            if flipxy
                println(io,"v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io,"v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            t = mesh.botTriangles[row, col]
            if reverse
                println(io, "f $(t.pt3) $(t.pt2) $(t.pt1)")
            else
                println(io, "f $(t.pt1) $(t.pt2) $(t.pt3)")
            end
        end

        println(io, "dims $(mesh.width) $(mesh.height)")
    end
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function stl_2_mesh(filename)
    lines = readlines(filename)

    vertexLines = [l for l in lines if startswith(l, "v")]
    nodeList = Vector{Point3D}(undef, size(vertexLines))

    count = 0
    for line in vertexLines
        count += 1
        elements = split(line, " ")
        x = parse(Float64, elements[2])
        y = parse(Float64, elements[3])
        z = parse(Float64, elements[4]) * 10
        pt = Point3D(x, y, z, 0, 0)
        nodeList[count] = pt
    end

    faceLines = [l for l in lines if startswith(l, "f")]
    triangles = Vector{Triangle}(undef, size(faceLines))
    for line in faceLines
        elements = split(line, " ")
        triangles[line] = Triangle(
            parse(Int64, elements[2]),
            parse(Int64, elements[3]),
            parse(Int64, elements[4]),
        )
    end

    dimsLines = [l for l in lines if startswith(l, "dims")]
    elements = split(dimsLines[1], " ")

    return Mesh(nodeList, triangles, parse(Int64, elements[2]), parse(Int64, elements[3]))
end


"""
$(SIGNATURES)
"""
function ∇(ϕ::Matrix{Float64})
    width, height = size(f)

    ∇ϕᵤ = zeros(Float64, width, height)   # the right edge will be filled with zeros
    ∇ϕᵥ = zeros(Float64, width, height)   # the buttom edge will be filled with zeros

    for x = 1:width-1, y = 1:height-1
        ∇ϕᵤ[x, y] = ϕ[x+1, y] - ϕ[x, y]
        ∇ϕᵥ[x, y] = ϕ[x, y+1] - ϕ[x, y]
    end

    return ∇ϕᵤ, ∇ϕᵥ
end


"""
$(SIGNATURES)
"""
function getPixelArea(mesh::Mesh)
    # A Mesh is a grid of 3D points. The X and Y coordinates are not necessarily aligned or square
    # The Z coordinate represents the value. brightness is just proportional to area.
    height, width = size(mesh.rectangles)

    pixelAreas = zeros(Float64, height, width)

    for row = 1:mesh.height, col = 1:width
        #=
        *------*
        |    / |
        |   /  |
        |  /   |
        | /    |
        *------*
        =#
        pixelAreas[width, col] =
            triangle_area(mesh.topTriangles[height, width]) +
            triangle_area(mesh.botTriangles[height, width])
    end

    return topPixelAreas
end


"""
$(SIGNATURES)
"""
function relax!(height_map::Matrix{Float64}, D::Matrix{Float64})
    # This function implements successive over relaxation for a matrix and its associated error matrix
    # There is a hardcoded assumption of Neumann boundary conditions--that the derivative across the
    # boundary must be zero in all cases. See:
    # https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

    width, height = size(height_map)
    val_average = sum(height_map) / (width * height)

    # Embed matrix within a larger matrix for better vectorization
    embedding = zeros(Float64, height + 2, width + 2)

    embedding[:, :] .= val_average
    embedding[2:height+1, 2:height+1] .= height_map[:, :]

    val_up =    embedding[1:height,   2:width+1]
    val_down =  embedding[3:height+2, 2:width+1]
    val_left =  embedding[2:height+1, 1:width]
    val_right = embedding[2:height+1, 3:width+2]

    delta = (val_up + val_down + val_left + val_right - D) ./ 4.0
    delta = ω .* (delta - height_map)
    height_map += delta

    max_update = maximum(abs.(delta))


    # node.z = .25 * (node_up.z + node_down.z + node_left.z + node_right.z) # simple averaging
    # node.z += ω/4 * (node_up.z + node_down.z + node_left.z + node_right.z - 4 * node.z)

    # matrix[x, y] += ω/4 * (val_up + val_down + val_left + val_right - 4 * val - D[x, y])

    return max_update
end


"""
$(SIGNATURES)

This function will take a `grid_definition x grid_definition` matrix and returns a
`grid_definition x grid_definition` mesh.
"""
function matrix_to_mesh(height_map::Matrix{Float64})
    height, width = size(height_map)

    result_mesh = create_mesh(height, width)
    for row = 1:height, col = 1:width
        result_mesh.rectangles[row, col].z = height_map[row, col]

        result_mesh.topPixelAreas[row, col].pt1.z = height_map[row, col]
        result_mesh.topNodes[row, col].z = height_map[row, col]

        result_mesh.botPixelAreas[row, col].pt1.z = height_map[row, col]
        result_mesh.botNodes[row, col].z = height_map[row, col]
    end

    return result_mesh
end


"""
$(SIGNATURES)
"""
function marchMesh!(mesh::Mesh, ϕ::Matrix{Float64})
    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)

    height, imgHeight = size(mesh.rectangles)

    # For each point in the mesh we need to figure out its velocity
    for x = 1:height, y = 1:mesh.height
        # XY coordinates in the mesh ARE XY coordinates in the image. The mesh just needs an extra row and column
        # at the bottom right edge so that the triangles can be closed (a triangle per pixel)

        if x == height
            u = 0
        else
            u = (y == mesh.height ? ∇ϕᵤ[x, y-1] : ∇ϕᵤ[x, y])
        end

        if y == mesh.height
            v = 0
        else
            v = (x == height ? ∇ϕᵥ[x-1, y] : ∇ϕᵥ[x, y])
        end

        velocities[x, y] = Point3D(-u, -v, 0, 0, 0)
    end

    # Basically infinity time to shrink a triangle to zip.
    min_t = 10_000.0

    for i ∈ 1:length(mesh.topTriangles)
        triangle = mesh.topTriangles[i]

        p1 = mesh.topNodes[triangle.pt1]
        p2 = mesh.topNodes[triangle.pt2]
        p3 = mesh.topNodes[triangle.pt3]

        v1 = velocities[p1.ix, p1.iy]
        v2 = velocities[p2.ix, p2.iy]
        v3 = velocities[p3.ix, p3.iy]

        # Get the time, at that velocity, for the area of the triangle to be nil.
        t1, t2 = find_maximum_t(p1, p2, p3, v1, v2, v3)

        # We are only interested in positive times.
        if !ismissing(t1)
            if 0 < t1 < min_t
                min_t = t1
            end

            if 0 < t2 < min_t
                min_t = t2
            end
        end
    end

    println("Overall min_t:", min_t)

    # Modify the mesh triangles but ensuring that we are far from destroying any of them.
    δ = min_t / 2


    for point in mesh.topNodes
        v = velocities[point.ix, point.iy]
        point.x = v.x * δ + point.x
        point.y = v.y * δ + point.y
    end

    # saveObj(mesh, "gateau.obj")
end


"""
$(SIGNATURES)
"""
function quantifyLoss!(D, suffix, img)
    println("Loss:")
    println("\tMinimum: $(minimum(D))")
    println("\tMaximum: $(maximum(D))")

    blue = zeros(size(D)); blue[D.>0] = D[D.>0]
    red = zeros(size(D));  red[D.<0] = -D[D.<0]
    green = zeros(size(D))

    println(size(blue))
    println(size(red))
    println(size(green))

    rgbImg = RGB.(red, green, blue)'
    save("./examples/loss_$(suffix).png", map(clamp01nan, rgbImg))

    # println("Saving output image:")
    # println(typeof(img))
    # E = Gray.(D)
    # println(typeof(E))
    # outputImg = img - E
    # save("./examples/actual_$(suffix).png", outputImg)
end


"""
$(SIGNATURES)
"""
function oneIteration(meshy, img, suffix)
    # Remember meshy is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.
    LJ = getPixelArea(meshy)
    D = Float64.(LJ - img)

    # Save the loss image as a png
    println(minimum(D))
    println(maximum(D))
    quantifyLoss!(D, suffix, img)

    # ∇Lᵤ, ∇Lᵥ = ∇(D)
    # plotVAsQuiver(∇Lᵤ, ∇Lᵥ, stride=10, scale=10, max_length=200)
    # println("okay")
    # return

    # save("./examples/loss_$(suffix).png", colorview(Gray, D))
    # return
    width, height = size(img)

    ϕ = Matrix{Float64}(undef, width, height)

    for i = 1:Grid_Definition * 10
        max_update = relax!(ϕ, D)

        i % 500 == 0 && println(max_update)

        max_update < 0.00001 &&
            println("Convergence reached at step $(i) with max_update of $(max_update)")
        break
    end

    save_stl!(
        matrix_to_mesh(ϕ * 0.02),
        "./examples/phi_$(suffix).obj",
        reverse = false,
        flipxy = true,
    )
    # plotAsQuiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reversex=false, reversey=false)
    # saveObj(matrix_to_mesh(D * 10), "D_$(suffix).obj")

    # Now we need to march the x,y locations in our mesh according to this gradient!
    marchMesh!(meshy, ϕ)
    save_stl!(meshy, "./examples/mesh_$(suffix).obj", flipxy = true)
end


"""
$(SIGNATURES)
"""
function setHeights!(mesh, heights, heightScale = 1.0, heightOffset = 50)
    width, height = size(heights)

    for y = 1:height, x = 1:width
        mesh.nodeArray[x, y].z = heights[x, y] * heightScale + heightOffset
        if x == 100 && y == 100
            println(
                "Example heights: $(heights[x, y]) and $(heights[x, y] * heightScale) and $(heights[x, y] * heightScale + heightOffset)",
            )
        end
    end

    # get the side edge
    for y = 1:height
        mesh.nodeArray[width+1, y].z = mesh.nodeArray[width, y].z
    end

    # get the bottom edge including the pesky corner
    for x = 1:width+1
        mesh.nodeArray[x, height+1].z = mesh.nodeArray[x, height].z
    end

    # # get the pesky corner!
    # mesh.nodeArray[width + 1, height + 1].z = mesh.nodeArray[width, height].z
end


"""
$(SIGNATURES)
"""
function setHeights(mesh, heightmap)
    width = mesh.width
    height = mesh.height
    nodes = Vector{Point3D}(undef, size(mesh.nodes))
    nodeArray = Matrix{Point3D}(undef, 0, 0)
    triangles = Vector{Triangle}(undef, size(mesh.triangles))
    scale = 1.0

    count = 0
    w, h = size(heightmap)
    for y = 1:height, x = 1:width
        count += 1
        point = mesh.nodes[count]

        z = heightmap[x, y]

        nodes[count] = Point3D(point.x * scale, point.y * scale, z * scale, 0, 0)
    end

    return Mesh(nodes, nodeArray, triangles, width, height)
end


"""
$(SIGNATURES)
"""
function solidify(inputMesh, offset = 100)
    width = inputMesh.width
    height = inputMesh.height
    totalNodes = width * height * 2

    nodeList = Vector{Point3D}(undef, totalNodes)
    nodeArrayTop = Matrix{Point3D}(undef, width, height)
    nodeArrayBottom = Matrix{Point3D}(undef, width, height)

    # imagine a 4x4 image. 4 * 2 + 2 * 2 = 12
    numEdgeNodes = width * 2 + (height - 2) * 2

    numTrianglesBottom = numTrianglesTop = (width - 1) * (height - 1) * 2
    numTrianglesEdges = numEdgeNodes * 2

    totalTriangles = numTrianglesBottom + numTrianglesTop + numTrianglesEdges

    println(
        "Specs: $(width)  $(height)  $(totalNodes)  $(numEdgeNodes)  $(numTrianglesBottom) $(totalTriangles)",
    )

    # Build the bottom surface
    count = 0
    for y = 1:height, x = 1:width
        count += 1
        newPoint = Point3D(x, y, -offset, x, y)
        nodeList[count] = newPoint
        nodeArrayBottom[x, y] = newPoint
    end

    # Copy in the top surface
    for y = 1:height, x = 1:width
        count += 1
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
    end

    println("We now have $(count-1) valid nodes")

    # Build the triangles for the bottom surface
    triangles = Vector{Triangle}(undef, totalTriangles)
    count = 0
    for y = 1:(height-1), x = 1:(width-1)
        count += 1
        # here x and y establish the column of squares we're in
        index_ul = (y - 1) * width + x
        index_ur = index_ul + 1

        index_ll = y * width + x
        index_lr = index_ll + 1

        triangles[2*count-1] = Triangle(index_ul, index_ll, index_ur)
        triangles[2*count] = Triangle(index_lr, index_ur, index_ll)
    end

    println("We've filled up $(count-1) triangles")
    if 2*count != numTrianglesBottom
        println(
            "Hmm aren't count and triangles bottom equal? $(count) vs $(numTrianglesBottom)",
        )
    end

    # Build the triangles for the top surface
    for y = 1:(height-1), x = 1:(width-1)
        count += 1

        # here x and y establish the column of squares we're in
        index_ul = (y - 1) * width + x + totalNodes / 2
        index_ur = index_ul + 1

        index_ll = y * width + x + totalNodes / 2
        index_lr = index_ll + 1

        triangles[2*count-1] = Triangle(index_ul, index_ur, index_ll)
        triangles[2*count] = Triangle(index_lr, index_ll, index_ur)
    end

    println("We've filled up $(count-1) triangles")

    # Build the triangles to close the mesh
    x = 1
    for y = 1:(height-1)
        count += 1

        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = y * width + x
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ul, ur)
        triangles[2*count] = Triangle(ur, lr, ll)
    end

    x = width
    for y = 1:(height-1)
        count += 1

        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = y * width + x
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ur, ul)
        triangles[2*count] = Triangle(ur, ll, lr)
    end

    y = 1
    for x = 2:width
        count += 1

        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = (y - 1) * width + (x - 1)
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ul, ur)
        triangles[2*count] = Triangle(ur, lr, ll)
    end

    y = height
    for x = 2:width
        count += 1

        ll = (y - 1) * width + x
        ul = ll + totalNodes / 2
        lr = (y - 1) * width + (x - 1)
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ur, ul)
        triangles[2*count] = Triangle(ur, ll, lr)
    end

    return Mesh(nodeList, nodeArrayBottom, triangles, width, height)
end


"""
$(SIGNATURES)
"""
function findSurface(mesh, image, f, imgWidth)
    width, height = size(image)

    # imgWidth = .1 # m
    # f = 1.0  # m
    H = f
    metersPerPixel = imgWidth / width
    println(metersPerPixel)

    Nx = zeros(Float64, width + 1, height + 1)
    Ny = zeros(Float64, width + 1, height + 1)

    for j = 1:height, i = 1:width
        node = mesh.nodeArray[i, j]
        dx = (node.ix - node.x) * metersPerPixel
        dy = (node.iy - node.y) * metersPerPixel

        little_h = node.z * metersPerPixel
        H_minus_h = H - little_h
        dz = H_minus_h


        # k = η * sqrt(dx * dx + dy * dy + H_minus_h * H_minus_h) - H_minus_h
        # Nx[i, j] = 1/k * dx
        # Ny[i, j] = 1/k * dy
        Ny[i, j] = tan(atan(dy / dz) / (n₁ - n₂))
        Nx[i, j] = tan(atan(dx / dz) / (n₁ - n₂))
    end


    # We need to find the divergence of the Vector field described by Nx and Ny
    divergence = zeros(Float64, width, height)

    for j = 1:height, i = 1:width
        δx = (Nx[i+1, j] - Nx[i, j])
        δy = (Ny[i, j+1] - Ny[i, j])
        divergence[i, j] = δx + δy

        i == 100 && j == 100 && println("div: $(divergence[i, j])")
    end
    println("Have all the divergences")

    h = zeros(Float64, width, height)
    max_update = 0
    for i = 1:Grid_Definition * 10 / 2
        max_update = relax!(h, divergence)

        i % 100 == 0 && println(max_update)
        max_update < 0.00001 &&
            println("Convergence reached at step $(i) with max_update of $(max_update)")
        break
    end
    # saveObj(matrix_to_mesh(h / 10), "./examples/heightmap.obj")
    return h, metersPerPixel
end


"""
$(SIGNATURES)
"""
function testSquareMesh!()
    mesh = create_mesh(100, 50)

    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodeArray[1, 1].x = 8
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])

    mesh.nodes[1].y += 12
    println(mesh.nodeArray[1, 1])
    println(mesh.nodes[1])
end


"""
$(SIGNATURES)
"""
function testSolidify!()
    println("Testing solidification")
    width = 100
    height = 100
    origMesh = create_mesh(width, height)

    for y = 1:height, x = 1:width
        x2 = (x - width / 2) / width
        y2 = (y - height / 2) / height
        value = x2 * x2 + y2 * y2
        origMesh.nodeArray[x, y].z = 15 - value * 25
    end

    save_stl!(origMesh, "./examples/testSolidify.obj")
    solidMesh = solidify(origMesh, 0)
    save_stl!(solidMesh, "./examples/testSolidify2.obj")
end


"""
$(SIGNATURES)
"""
function plotAsQuiver(
    g;
    stride = 4,
    scale = 300,
    max_length = 2,
    flipxy = false,
    reversey = false,
    reversex = false,
)

    h, w = size(g)
    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]

    for x = 1:stride:w, y = 1:stride:h
        reversex ? push!(xs, x) : push!(xs, -x)
        reversey ? push!(ys, -y) : push!(ys, y)

        p1 = g[y, x]
        u = (g[y, x+1] - g[y, x]) * scale
        v = (g[y+1, x] - g[y, x]) * scale

        u = -u

        reversey && (v = -v)
        reversex && (u = -u)

        # println(u, v)
        u >= 0 ? push!(us, min(u, max_length)) : push!(us, max(u, -max_length))
        v >= 0 ? push!(vs, min(v, max_length)) : push!(vs, max(v, -max_length))
    end

    q =
        flipxy ? quiver(ys, xs, quiver = (vs, us), aspect_ratio = :equal) :
        quiver(xs, ys, quiver = (us, vs), aspect_ratio = :equal)

    display(q)
    readline()
end


"""
$(SIGNATURES)
"""
function plotVAsQuiver(vx, vy; stride = 4, scale = 300, max_length = 2)
    h, w = size(vx)

    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]

    for x = 1:stride:w, y = 1:stride:h
        push!(xs, x)
        push!(ys, h - y)

        u = max(vx[x, y], 0.001)
        v = max(vy[x, y], 0.001)

        push!(us, u)
        push!(vs, v)
        # println(u, ": ", v)
    end

    # readline()
    q = quiver(xs, ys, quiver = (us, vs), aspect_ratio = :equal)
    display(q)
    readline()
end



"""
$(SIGNATURES)
"""
function engineer_caustics(img)
    img = Gray.(img)
    img2 = permutedims(img) * 1.0
    width, height = size(img2)

    # meshy is the same size as the image with an extra row/column to have coordinates to
    # cover each image pixel with a triangle.
    meshy = create_mesh(width + 1, height + 1)

    # We need to boost the brightness of the image so that its sum and the sum of the area are equal
    mesh_sum = width * height
    image_sum = sum(img2)
    boost_ratio = mesh_sum / image_sum

    # img3 is `grid_definition x grid_definition` and is normalised to the same (sort of) _energy_ as the
    # original image.
    img3 = img2 .* boost_ratio

    oneIteration(meshy, img3, "it1")
    oneIteration(meshy, img3, "it2")
    oneIteration(meshy, img3, "it3")

    # oneIteration(meshy, img3, "it4")
    # oneIteration(meshy, img3, "it5")
    # oneIteration(meshy, img3, "it6")

    h, metersPerPixel = findSurface(meshy, img3, 1.0, Artifact_Size)

    setHeights!(meshy, h, 1.0)
    # newMesh = setHeights(meshy, h)

    solidMesh = solidify(meshy)
    save_stl!(
        solidMesh,
        "./examples/original_image.obj",
        scale = Float64(1 / Grid_Definition * Artifact_Size),
        scalez = Float64(1 / Grid_Definition * Artifact_Size)
    )

    return meshy, img3
end


"""
$(SIGNATURES)
"""
function main()
    @assert size(ARGS) == (1,) "Intented usage is: julia create_mesh.jl image.png"

    img = Images.load(ARGS[1])
    return engineer_caustics(img)
end


"""
$(SIGNATURES)
"""
main()
