"""
$(SIGNATURES)

Given 3 points and their velocities, calculate the time `t` required to bring the area of that triangle to zero
"""
function find_maximum_t(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D)
    # The area of a triangle is 1/2 * [ Ax (By - Cy) + Bx (Cy - Ay) + Cx (Ay - By)]
    # where each point of the triangle is where it will be after time t
    # i.e. a point goes from P to P+tV where V is the velocity of that point.
    # Notation is here B -> B + t_vB

    # To make the calculation simpler, everything is translated so that A is at the
    # origin of the plane and its velocity is nil.
    Bx = p2.x - p1.x
    By = p2.y - p1.y
    Cx = p3.x - p1.x
    Cy = p3.y - p1.y

    t_vBx = p2.vx - p1.vx
    t_vBy = p2.vy - p1.vy
    t_vCx = p3.vx - p1.vx
    t_vCy = p3.vy - p1.vy

    # After this, given that Ax = Ay = t_vAx = t_vAy = 0, the area is nil iff
    # (Bx + t_vBx) (Cy + t_vCy )- (Cx + t_vCx) (By + t_vBy) = 0.
    # After expansion and reshuffling to have a quadratic equation where t
    # is the variable, the coefficients of that equation are:
    a = t_vBy * t_vCx - t_vCy * t_vBx
    b = By * t_vCx + Cx * t_vBy - (Bx * t_vCy + Cy * t_vBx)
    c = Cx * By - Bx * Cy

    # if a = 0, this is just a linear equation.
    if a == 0
        # Fingers crossed that b != 0 but extremely unlikely since B is a Float64.
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
    mesh::RectangleMesh,
    filename::String;
    scale = 1.0,
    scalez = 1.0,
    reverse = false,
    flipxy = false,
)

    height, width = size(mesh)

    open(filename, "w") do io
        for row ∈ 1:height, col ∈ 1:width

            # Top triangle
            top_triangle = top_triangle(mesh, row, col)
            n = centroid(top_triangle)

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end


            # Bottom triangle
            bot_triangle = bottom_triangle(mesh, row, col)
            n = centroid(bot_triangle)

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end

            if flipxy
                println(io, "v $(n.y * scale) $(n.x * scale) $(n.z * scalez)")
            else
                println(io, "v $(n.x * scale) $(n.y * scale) $(n.z * scalez)")
            end


        end

        # CHECK what dims exactly represents. Number of triangles?
        println(io, "dims $(2*mesh.width) $(2*mesh.height)")
    end
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function stl_2_mesh(filename)
    lines = readlines(filename)

    vertexLines = [l for l in lines if startswith(l, "v")]
    nodeList = Vector{Vertex3D}(undef, size(vertexLines))

    count = 0
    for line in vertexLines
        count += 1
        elements = split(line, " ")
        x = parse(Float64, elements[2])
        y = parse(Float64, elements[3])
        z = parse(Float64, elements[4]) * 10
        pt = Vertex3D(x, y, z, 0, 0)
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

    return RectangleMesh(
        nodeList,
        triangles,
        parse(Int64, elements[2]),
        parse(Int64, elements[3]),
    )
end


"""
$(SIGNATURES)
"""
function ∇(ϕ::Matrix{Float64})
    width, height = size(f)

    ∇ϕᵤ = zeros(Float64, width, height)   # the right edge will be filled with zeros
    ∇ϕᵥ = zeros(Float64, width, height)   # the buttom edge will be filled with zeros

    for row = 1:height, col = 1:width
        ∇ϕᵤ[row, col] = ϕ[row+1, col] - ϕ[row, col]
        ∇ϕᵥ[row, col] = ϕ[row, col+1] - ϕ[row, col]
    end

    return ∇ϕᵤ, ∇ϕᵥ
end


"""
$(SIGNATURES)
"""
function get_pixel_area(mesh::RectangleMesh)
    # A Mesh is a grid of 3D points. The X and Y coordinates are not necessarily aligned or square
    # The Z coordinate represents the value. brightness is just proportional to area.
    height, width = size(mesh)

    pixelAreas = zeros(Float64, height, width)

    for row = 1:height, col = 1:width
        #
        # *------*
        # | TOP/ |
        # |   /  |
        # |  /   |
        # | /BOT |
        # *------*
        #
        pixelAreas[row, col] =
            triangle_area(top_triangles(mesh, row, col)) +
            triangle_area(bottom_triangle(mesh, row, col))
    end

    return pixelAreas
end


"""
$(SIGNATURES)

This function implements successive over relaxation for a matrix and its associated error matrix
There is a hardcoded assumption of Neumann boundary conditions--that the derivative across the
boundary must be zero in all cases. See:
https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

"""
function relax!(height_map::Matrix{Float64}, divergence::Matrix{Float64})
    height, width = size(height_map)
    n_pixels = height * width
    val_average = sum(height_map) / n_pixels

    # Embed matrix within a larger matrix for better vectorization and avoid duplicated code
    container = zeros(Float64, height + 2, width + 2)

    container[:, :] .= val_average
    container[2:height+1, 2:height+1] .= height_map[:, :]

    val_up = container[1:height, 2:width+1]
    val_down = container[3:height+2, 2:width+1]
    val_left = container[2:height+1, 1:width]
    val_right = embcontaineredding[2:height+1, 3:width+2]

    delta = (val_up + val_down + val_left + val_right - divergence) ./ 4.0
    height_map = height_map + ω .* (delta - height_map)

    return maximum(abs.(delta))
end


"""
$(SIGNATURES)

This function will take a `grid_definition x grid_definition` matrix and returns a
`grid_definition x grid_definition` mesh.

"""
function matrix_to_mesh(height_map::Matrix{Float64})
    height, width = size(height_map)

    mesh = RectangleMesh(height, width)

    for row ∈ 1:height, col ∈ 1:width
        mesh.vertexList[mesh.to_index(row, col)].z = height_map[row, col]
    end

    # The borders height is forced at 0.
    for r ∈ 1:height+1
        mesh.vertexList[mesh.to_index(r, 1)].z = 0.0
        mesh.vertexList[mesh.to_index(r, width + 1)].z = 0.0
    end

    for c ∈ 1:width+1
        mesh.vertexList[mesh.to_index(1, c)].z = 0.0
        mesh.vertexList[mesh.to_index(height + 1, c)].z = 0.0
    end

    return mesh
end


"""
$(SIGNATURES)
"""
function march_mesh!(mesh::RectangleMesh, ϕ::Matrix{Float64})
    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)

    height, width = size(mesh)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.
    for row ∈ 1:height, col ∈ 1:width
        Vx = ∇ϕᵤ[x, y]
        Vy = ∇ϕᵥ[x, y]

        mesh.vertexList[mesh.to_index(row, col)].vx = -Vx
        mesh.vertexList[mesh.to_index(row, col)].vy = -Vy
    end

    # The velocity matrix was initialized with zeros everywhere. We nevertheless ovewrite them just in case...
    mesh.rectangles[:, 1].vx .= 0.0
    mesh.rectangles[:, end].vy = 0.0

    mesh.rectangles[1, :].vx .= 0.0
    mesh.rectangles[end, :].vy = 0.0


    # Basically infinity time to shrink a triangle to zip.
    min_t = 10_000.0

    for row ∈ 1:height, col ∈ 1:width
        top_tri = top_triangles(mesh, row, col)
        p1 = top_tri.pt1
        p2 = top_tri.pt2
        p3 = top_tri.pt3

        # Get the time, at that velocity, for the area of the triangle to be nil.
        # We are only interested in positive times.
        t1, t2 = find_maximum_t(p1, p2, p3)
        for t ∈ [t1, t2]
            if (!ismissing(t)) && (0 < t < min_t)
                min_t = t
            end
        end


        bottom_tri = bottom_triangle(mesh, row, col)
        p1 = bottom_tri.pt1
        p2 = bottom_tri.pt2
        p3 = bottom_tri.pt3

        # Get the time, at that velocity, for the area of the triangle to be nil.
        # We are only interested in positive times.
        t1, t2 = find_maximum_t(p1, p2, p3)
        for t ∈ [t1, t2]
            if (!ismissing(t)) && (0 < t < min_t)
                min_t = t
            end
        end

    end

    println("Overall min_t:", min_t)

    # Modify the mesh triangles but ensuring that we are far from destroying any of them with a nil area.
    δ = min_t / 2

    for row ∈ 1:height, col ∈ 1:width
        Vx = mesh.rectangles[row, col].vx
        Vy = mesh.rectangles[row, col].vy

        mesh.rectangles[row, col].x += δ * Vx
        mesh.rectangles[row, col].y += δ * Vy
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

    blue = zeros(size(D))
    blue[D.>0] = D[D.>0]
    red = zeros(size(D))
    red[D.<0] = -D[D.<0]
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
function single_iteration!(mesh, image, suffix)
    # Remember meshy is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.
    LJ = get_pixel_area(mesh)
    D = Float64.(LJ - image)

    # Save the loss image as a png
    println(minimum(D))
    println(maximum(D))
    quantifyLoss!(D, suffix, image)

    width, height = size(image)

    ϕ = Matrix{Float64}(undef, width, height)

    for i = 1:N_Iterations_Convergence
        max_update = relax!(ϕ, D)

        i % 500 == 0 && println(max_update)

        max_update < 0.00001 && (
            println("Convergence reached at step $(i) with max_update of $(max_update)"); break
        )
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
    march_mesh!(mesh, ϕ)
    save_stl!(mesh, "./examples/mesh_$(suffix).obj", flipxy = true)
end


"""
$(SIGNATURES)
"""
function set_heights!(
    mesh,
    height_map,
    height_scale = 1.0,
    height_offset = Top_Offset / Meters_Per_Pixel,
)
    width, height = size(height_map)

    for row ∈ 1:height, col ∈ 1:width
        mesh.rectangles[row, col].z = height_map[x, y] * height_scale + height_offset

        if x == 100 && y == 100
            println(
                "Example heights: $(height_map[x, y]) and $(height_map[x, y] * heightScale) and $(height_map[x, y] * heightScale + heightOffset)",
            )
        end
    end

    # Forces the borders at Height_Offset. They will never change because nil velocity.
    mesh.rectangles[1, :] .= height_offset
    mesh.rectangles[end, :] .= height_offset
    mesh.rectangles[:, 1] .= height_offset
    mesh.rectangles[:, end] .= height_offset
end


"""
$(SIGNATURES)
"""
function set_heights(height_map, height_scale = 1.0)
    height, width = size(height_map)
    mesh = RectangleMesh(height, width)

    return set_heights!(mesh, height_map, height_scale = 1.0)
end


"""
$(SIGNATURES)
"""
function create_solid(mesh, bottom_offset = Bottom_Offset / Meters_Per_Pixel)
    height, width = size(mesh)

    ###
    ### Top mesh
    ###
    # List top mesh top triangles
    list_triangles_top_top = [top_triangles(mesh, r, c) for r ∈ 1:height, c ∈ 1:width]

    # List top mesh bottom triangles
    list_triangles_top_bottom = [bottom_triangle(mesh, r, c) for r ∈ 1:height, c ∈ 1:width]


    ###
    ### Botttom mesh
    ###
    # Build the bottom surface
    bottom_mesh = RectangleMesh(height, width)
    for row = 1:height+1, col = 1:width+1
        bottom_mesh.rectangles[row, col] = Vertex3D(row, col, -offset, 0.0, 0.0)
    end

    # List top mesh top triangles
    list_triangles_bot_top =
        [top_triangles(bottom_mesh, r, c) for r ∈ 1:height, c ∈ 1:width]

    # List top mesh bottom triangles
    list_triangles_bot_bottom =
        [bottom_triangle(bottom_mesh, r, c) for r ∈ 1:height, c ∈ 1:width]



    ###
    ### Side meshes
    ###
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

    return RectangleMesh(nodeList, nodeArrayBottom, triangles, width, height)
end


"""
$(SIGNATURES)
"""
function find_surface(mesh, image, f = Focal_Length, picture_width = Picture_Side)
    height, width = size(mesh)

    H = f
    metersPerPixel = picture_width / width
    println(metersPerPixel)

    # h(eight) indexes rows, w(idth) indexes columns
    Nx = zeros(Float64, height + 1, width + 1)
    Ny = zeros(Float64, height + 1, width + 1)

    for row ∈ 1:height, col ∈ 1:width
        node = mesh.nodeArray[row, col]
        dx = (node.ix - node.x) * metersPerPixel
        dy = (node.iy - node.y) * metersPerPixel

        little_h = node.z * metersPerPixel
        H_minus_h = H - little_h
        dz = H_minus_h

        # k = η * sqrt(dx * dx + dy * dy + H_minus_h * H_minus_h) - H_minus_h
        Nx[row, col] = tan(atan(dx / dz) / (n₁ - n₂))
        Ny[row, col] = tan(atan(dy / dz) / (n₁ - n₂))
    end


    # We need to find the divergence of the Vector field described by Nx and Ny
    divergence = zeros(Float64, width, height)

    for row ∈ 1:height, col ∈ 1:width
        δx = Nx[row+1, col] - Nx[row, col]
        δy = Ny[row, col+1] - Ny[row, col]
        divergence[row, col] = δx + δy

        (row % 100 == 0) &&
            (col % 100 == 0) &&
            println("div[$(row), $(col)]: $(divergence[row, col])")
    end
    println("Have all the divergences")

    height_map = zeros(Float64, height, width)
    max_update = 0.0
    for iter = 1:N_Iterations_Convergence
        max_update = relax!(height_map, divergence)

        iter % 100 == 0 && println(max_update)
        if max_update < 0.00001
            println("Convergence reached at step $(iter) with max_update of $(max_update)")
            break
        end
    end

    # saveObj(matrix_to_mesh(h / 10), "./examples/heightmap.obj")
    return height_map, metersPerPixel
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
    mesh = RectangleMesh(width + 1, height + 1)

    # We need to boost the brightness of the image so that its sum and the sum of the area are equal
    mesh_sum = width * height
    image_sum = sum(img2)
    boost_ratio = mesh_sum / image_sum

    # img3 is `grid_definition x grid_definition` and is normalised to the same (sort of) _energy_ as the
    # original image.
    img3 = img2 .* boost_ratio

    single_iteration!(mesh, img3, "it1")
    single_iteration!(mesh, img3, "it2")
    single_iteration!(mesh, img3, "it3")

    # single_teration(meshy, img3, "it4")
    # single_teration(meshy, img3, "it5")
    # single_teration(meshy, img3, "it6")


    h, metersPerPixel = find_surface(mesh, img3, 1.0, Artifact_Size)

    set_heights!(mesh, h)

    solidMesh = create_solid(mesh)
    save_stl!(
        solidMesh,
        "./examples/original_image.obj",
        scale = Float64(1 / Grid_Definition * Artifact_Size),
        scalez = Float64(1 / Grid_Definition * Artifact_Size),
    )

    return mesh, img3
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
