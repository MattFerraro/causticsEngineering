"""
$(SIGNATURES)
"""
function engineer_caustics(source_image)
    imageBW = Float64.(Gray.(source_image))

    width, height = size(imageBW)
    println("Image size: $((height, width))")

    # mesh is the same size as the image with an extra row/column to have coordinates to
    # cover each image pixel with a triangle.
    mesh = FaceMesh(height, width)
    print("Mesh creation: ")
    println("$(mesh.topleft.vx[15, 5]) , $(mesh.topleft.vy[15, 5])")

    # We need to boost the brightness of the image so that its sum and the sum of the area are equal
    average_intensity = sum(imageBW) / (width * height)

    # imageBW is `grid_definition x grid_definition` and is normalised to the same (sort of) _energy_ as the
    # original image.
    imageBW = imageBW ./ average_intensity

    for i ∈ 1:2
        println("\nSTARTING ITERATION $(i) ---")
        max_update = single_iteration!(mesh, imageBW, "it$(i)")
        print("Single iteration: ")
        println("$(mesh.topleft.vx[15, 5]) , $(mesh.topleft.vy[15, 5])")

        println("---------- ITERATION $(i): $(max_update)")
    end

    height_map = find_surface(mesh, imageBW; f = 1.0, picture_width = Caustics_Side)
    print("Find surface: ")
    println("$(mesh.topleft.vx[15, 5]) , $(mesh.topleft.vy[15, 5])")

    set_heights!(mesh, height_map)
    print("Set heights!: ")
    println("$(mesh.topleft.vx[15, 5]) , $(mesh.topleft.vy[15, 5])")


    # solidMesh = create_solid(mesh)
    # save_stl!(
    #     solidMesh,
    #     "./examples/original_image.obj",
    #     scale = Float64(1 / Grid_Definition * Artifact_Size),
    #     scalez = Float64(1 / Grid_Definition * Artifact_Size),
    # )

    return mesh, imageBW
end



"""
$(SIGNATURES)
"""
function single_iteration!(mesh, image, suffix)
    # Remember mesh is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.

    # The idea is that:
    # - any pixel on the caustic projection receives light from a given 'rectangle' on the lens.
    # - That rectangle is made of 2 triangles.

    # DOES NOT WORK. D IS ALWAYS IDENTICAL
    LJ = get_area_pixels(mesh)
    D = Float64.(LJ - image)

    # Save the loss image as a png
    plot_loss!(D, suffix, image)

    # ϕ is the _heightmap_
    # ϕ = zeros(Float64, size(image))
    height_map = rand(Float64, size(image)) ./ 1_000 .- 0.5 / 1_000

    max_update = 100.0
    interaction_count = 0
    while true
        interaction_count += 1
        interaction_count % 500 == 0 && println("Converging for intensity: $(max_update)")

        old_max_update = max_update
        max_update = relax!(height_map, D)
        if abs(max_update) < 1e-6 ||
           abs((max_update - old_max_update) / old_max_update) < 0.01
            println(
                "Convergence stopped at step $(interaction_count) with max_update of $(max_update)",
            )
            break
        end
    end

    save_stl!(
        matrix_to_mesh(height_map * 0.02),
        "./examples/phi_$(suffix).obj",
        reverse = false,
        flipxy = true,
    )
    # plotAsQuiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reversex=false, reversey=false)
    # saveObj(matrix_to_mesh(D * 10), "D_$(suffix).obj")

    # Now we need to march the x,y locations in our mesh according to this gradient!
    march_mesh!(mesh, height_map)
    save_stl!(mesh, "./examples/mesh_$(suffix).obj", flipxy = true)

    return max_update
end


"""
$(SIGNATURES)
"""
function find_surface(
    mesh::FaceMesh,
    image;
    f = Focal_Length,
    picture_width = Caustics_Side,
)

    # height indexes rows, width indexes columns
    height, width = size(mesh)

    # Coordinates difference
    # C; Coordinates on the caustics. Simple rectangular values in pixels.
    # Coordinates on the lens face coming from topleft field.
    # x = rows; y = columns
    dx = repeat(Float64.(1:height+1), 1, width + 1) - mesh.topleft.x
    dy = repeat(Float64.(1:width+1)', height + 1, 1) - mesh.topleft.y
    height_map = height_map + ω .* (target_map - height_map)

    # true_H is the real distance beteen the surface of the lens and the projection
    # H is the initial distance from the surface of the carved face to the projection screen.
    # H is constant and does not account for the change in heights due to carving.
    # The focal length is in meters. H is in Pixels.
    H = f / Meters_Per_Pixel

    true_H = similar(mesh.topleft.z)
    @. true_H[:] = -mesh.topleft.z[:] + H

    # Normals.
    Nx = zeros(Float64, height + 1, width + 1)
    Ny = zeros(Float64, height + 1, width + 1)
    @. Nx[:] = tan(atan(dx[:] / true_H[:]) / (n₁ - n₂))
    @. Ny[:] = tan(atan(dy[:] / true_H[:]) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    divergence = zeros(Float64, width, height)
    δr = zeros(Float64, width, height)
    δc = zeros(Float64, width, height)

    @. δr[1:end, 1:end] = Nx[2:end, 1:end-1] - Nx[1:end-1, 1:end-1]
    @. δc[1:end, 1:end] = Nx[1:end-1, 2:end] - Ny[1:end-1, 1:end-1]
    divergence = δr + δc

    println("Have all the divergences")

    height_map = zeros(Float64, height, width)
    max_update = 1000.0
    iteration_count = 0
    while true
        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for divergence: $(max_update)")

        old_max_update = max_update
        max_update = relax!(height_map, divergence)
        if abs(max_update) < 1e-6 ||
           abs((max_update - old_max_update) / old_max_update) < 0.01
            println(
                "Convergence stopped improving at step $(iteration_count) with max_update of $(max_update)",
            )
            break
        end
    end

    # saveObj(matrix_to_mesh(h / 10), "./examples/heightmap.obj")
    return height_map
end


"""
$(SIGNATURES)

TODO: Think about the optimal offset for convergence.
"""
function set_heights!(
    mesh,
    height_map,
    height_scale = 1.0,
    height_offset = 0.0 * Top_Offset / Meters_Per_Pixel,
)
    width, height = size(height_map)

    @. mesh.topleft.z[1:height, 1:width] =
        height_map[1:height, 1:width] * height_scale + height_offset

    # Forces the borders at Height_Offset. They will never change because nil velocity.
    fill_borders!(mesh.topleft.z, height_offset)
end


"""
$(SIGNATURES)
"""
function set_heights(height_map, height_scale = 1.0)
    height, width = size(height_map)
    mesh = FaceMesh(height, width)

    return set_heights!(mesh, height_map, height_scale = 1.0)
end


"""
$(SIGNATURES)

A Mesh is a collection of triangles. The brightness flowing through a given triangle is just proportional to its
area in the x, y plane. z is ignored.

The function returns a matrix with the quantity of light coming from each 'rectangle'  around a pixel. That 'rectangle'
has been shifted and flexed around.
"""
function get_area_pixels(mesh::FaceMesh)
    height, width = size(mesh)

    pixel_areas = zeros(Float64, size(mesh))

    for ci ∈ CartesianIndices(pixel_areas)
        pixel_areas[ci] =
            area(top_triangle3D(mesh, ci)...) + area(bot_triangle3D(mesh, ci)...)
    end

    return pixel_areas
end


"""
$(SIGNATURES)
"""
function quantifyLoss!(D, suffix, img)
    println("Loss:")
    println("\tMinimum loss: $(minimum(D))")
    println("\tMaximum loss: $(maximum(D))")

    blue = zeros(size(D))
    blue[D.>0] = D[D.>0]
    red = zeros(size(D))
    red[D.<0] = -D[D.<0]
    green = zeros(size(D))

    # println(size(blue))
    # println(size(red))
    # println(size(green))

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

    container .= val_average
    container[2:height+1, 2:height+1] .= height_map[:, :]

    val_up = container[1:height, 2:width+1]
    val_down = container[3:height+2, 2:width+1]
    val_left = container[2:height+1, 1:width]
    val_right = container[2:height+1, 3:width+2]

    # Target position. The target is the current height map smoothed by averaging to which the
    # divergence is added.
    target_map = (val_up + val_down + val_left + val_right) ./ 4.0
    target_map += divergence

    # Let the heightmap converge towards the target at a slow rate.
    height_map = height_map + ω .* (target_map - height_map)

    max_update = maximum(abs.(height_map))
    println("Minimum / Maximum of the height map after relaxation.")
    println("\tMax update = $( maximum(height_map) )")
    println("\tMin update = $( minimum(height_map) )")

    return max_update
end


"""
$(SIGNATURES)

Given 3 points and their velocities, calculate the time `t` required to bring the area of that triangle to zero
"""
function find_maximum_t(p1::Vertex3D, p2::Vertex3D, p3::Vertex3D)
    # Three points A, B and C, with coordinates (x, y)
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
    a = t_vCy * t_vBx - t_vBy * t_vCx
    b = -By * t_vCx - Cx * t_vBy + Bx * t_vCy + Cy * t_vBx
    c = Bx * Cy - Cx * By

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
            # Very unlikely with Float64.
            return missing, missing
        end
    end
end

find_maximum_t(p::Tuple{Vertex3D,Vertex3D,Vertex3D}) = find_maximum_t(p[1], p[2], p[3])


"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh!(mesh::FaceMesh, ϕ::Matrix{Float64})

    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)
    println("Min/Max of values of ϕ: $(minimum(ϕ)) / $(maximum(ϕ))")
    println("Sum of values of ∇ϕᵤ:   $( sum(∇ϕᵤ) )")
    println("Sum of values of ∇ϕᵥ:   $( sum(∇ϕᵥ) )")

    height, width = size(ϕ)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.
    # Warning: The indices are necessary because topleft and ∇ϕ are of different sizes
    mesh.topleft.vx[1:height, 1:width] .= -∇ϕᵤ[1:height, 1:width]
    mesh.topleft.vy[1:height, 1:width] .= -∇ϕᵥ[1:height, 1:width]

    # The velocity matrix was initialized with zeros everywhere. We nevertheless ovewrite them just in case...
    fill_borders!(mesh.topleft.vx, 0.0)
    fill_borders!(mesh.topleft.vy, 0.0)

    # Basically infinity time to shrink a triangle to zip.
    min_positive_t = Inf
    max_negative_t = -Inf

    for row ∈ 1:height-1, col ∈ 1:width-1
        top_tri = top_triangle3D(mesh, row, col)

        # Get the time, at that velocity, for the area of the triangle to be nil.
        # We are only interested in positive times.
        t1, t2 = find_maximum_t(top_tri)
        for t ∈ [t1, t2]
            if (!ismissing(t)) && (0 < t < min_positive_t)
                min_positive_t = t
            end
            if (!ismissing(t)) && (max_negative_t < t < 0)
                max_negative_t = t
            end
        end

        bot_tri = bot_triangle3D(mesh, row, col)

        # Get the time, at that velocity, for the area of the triangle to be nil.
        # We are only interested in positive times.
        t1, t2 = find_maximum_t(bot_tri)
        for t ∈ [t1, t2]
            if (!ismissing(t)) && (0 < t < min_positive_t)
                min_positive_t = t
            end
            if (!ismissing(t)) && (max_negative_t < t < 0)
                max_negative_t = t
            end
        end
    end

    # Modify the mesh triangles but ensuring that we are far from destroying any of them with a nil area.
    # Use the maximum change possible (negative or positive)
    δ = (abs(max_negative_t) < min_positive_t) ? min_positive_t / 2.0 : max_negative_t / 2.0

    println("Overall maximum variations:")
    println("\tOverall min_positive_t: $(min_positive_t)")
    println("\tOverall max_negative_t: $(max_negative_t)")
    println("\tUsing δ: $(δ)")

    mesh.topleft.x[1:height, 1:width] .-= δ .* ∇ϕᵤ[1:height, 1:width]
    mesh.topleft.y[1:height, 1:width] .-= δ .* ∇ϕᵥ[1:height, 1:width]

    # saveObj(mesh, "gateau.obj")
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function create_solid(
    mesh::FaceMesh;
    bottom_offset = Bottom_Offset / Meters_Per_Pixel,
    top_offset = Top_Offset / Meters_Per_Pixel,
)
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
    # Build the bottom surface which is prepopulated by its constructor. However, its height is incorrect.

    # Coordinates on the caustics
    rows = repeat(Float64.(1:height), 1, width)
    cols = repeat(Float64.(1:width)', height, 1)

    bottom_mesh = FaceMesh(height, width)
    bottom_mesh.topleft.x[:] .= rows[:]
    bottom_mesh.topleft.y[:] .= cols[:]
    bottom_mesh.topleft.z[:] .= -Bottom_Offset
    bottom_mesh.topleft.vx[:] .= 0.0
    bottom_mesh.topleft.vx[:] .= 0.0


    ###
    ### Side meshes
    ###
    # Build triangles to create side meshes and close the mesh
    # The mesh is made of a single top/bottom triangles joining the bottom face
    # and the top face.

    # (Looking from above), left side.

    list_triangles = []
    col = 1
    for row = 1:height
        count += 1

        ll = (row - 1) * width + col
        ul = ll + totalNodes / 2
        lr = row * width + col
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ul, ur)
        triangles[2*count] = Triangle(ur, lr, ll)
    end

    mesh_right = FaceMesh(height, 1)
    for row = 1:(height-1)
        count += 1

        ll = (row - 1) * width + col
        ul = ll + totalNodes / 2
        lr = row * width + col
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ur, ul)
        triangles[2*count] = Triangle(ur, ll, lr)
    end

    mesh_up = FaceMesh(1, width)
    row = 1
    for col = 2:width
        count += 1

        ll = (row - 1) * width + col
        ul = ll + totalNodes / 2
        lr = (row - 1) * width + (col - 1)
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ul, ur)
        triangles[2*count] = Triangle(ur, lr, ll)
    end

    mesh_down = FaceMesh(1, width)
    for col = 2:width
        count += 1

        ll = (row - 1) * width + col
        ul = ll + totalNodes / 2
        lr = (row - 1) * width + (col - 1)
        ur = lr + totalNodes / 2
        triangles[2*count-1] = Triangle(ll, ur, ul)
        triangles[2*count] = Triangle(ur, ll, lr)
    end

    return RectangleMesh(nodeList, nodeArrarowBottom, triangles, width, height)
end




"""
$(SIGNATURES)
"""
function ∇(ϕ::Matrix{Float64})
    height, width = size(ϕ)

    ∇ϕᵤ = zeros(Float64, height, width)   # the right edge will be filled with zeros
    ∇ϕᵥ = zeros(Float64, height, width)   # the bottom edge will be filled with zeros

    for row ∈ 1:height-1, col ∈ 1:width-1
        ∇ϕᵤ[row, col] = ϕ[row+1, col] - ϕ[row, col]
        ∇ϕᵥ[row, col] = ϕ[row, col+1] - ϕ[row, col]
    end

    ∇ϕᵤ[1:end, end] .= 0.0
    ∇ϕᵤ[end, 1:end] .= 0.0

    ∇ϕᵥ[1:end, end] .= 0.0
    ∇ϕᵥ[end, 1:end] .= 0.0

    return ∇ϕᵤ, ∇ϕᵥ
end


"""
$(SIGNATURES)

This function will take a `grid_definition x grid_definition` matrix and returns a
`grid_definition x grid_definition` mesh.

"""
function matrix_to_mesh(height_map::Matrix{Float64})
    height, width = size(height_map)

    mesh = FaceMesh(height, width)

    # 1 more topleft than pixels! Therefore compiler needs to specify exact indices.
    mesh.topleft.z[1:height, 1:width] .= height_map[1:height, 1:width]

    # The borders' height is forced at 0.
    fill_borders!(mesh.topleft.z, 0.0)

    return mesh
end
