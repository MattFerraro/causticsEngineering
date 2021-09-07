"""
$(SIGNATURES)
"""
function engineer_caustics(source_image)
    imageBW = Float64.(Gray.(source_image))

    height, width = size(imageBW)
    println("Image size: $((height, width))")

    # mesh is the same size as the image with an extra row/column to have coordinates to
    # cover each image pixel with a triangle.
    mesh = FaceMesh(height, width)
    println("Mesh creation: ")
    r_check = 1
    c_check = 1
    println(
        "\tCoordinates: $(mesh.topleft.r[r_check, c_check]) , $(mesh.topleft.c[r_check, c_check]) , $(mesh.topleft.h[r_check, c_check])",
    )
    println(
        "\tVelocities: $(mesh.topleft.vr[r_check, c_check]) , $(mesh.topleft.vc[r_check, c_check])",
    )

    # The energy going through the lens is equal to the amount of energy on the caustics
    total_energy_lens = height * width
    total_energy_caustics = sum(imageBW)
    correction_ratio = (width * height) / sum(imageBW)

    # imageBW is `grid_definition x grid_definition` and is normalised to the same (sort of) _energy_ as the
    # original image.
    imageBW = imageBW .* correction_ratio

    for i ∈ 1:256
        println("\nSTARTING VERTICAL ITERATION $(i) ---")

        max_update = move_vertically!(mesh, imageBW, "it$(i)")
        print("Vertical move max update = $(max_update)")

        println(
            "\tCoordinates: $(mesh.topleft.r[r_check, c_check]) , $(mesh.topleft.c[r_check, c_check]) , $(mesh.topleft.h[r_check, c_check])",
        )
        println(
            "\tVelocities: $(mesh.topleft.vr[r_check, c_check]) , $(mesh.topleft.vc[r_check, c_check])",
        )

        println("---------- ITERATION $(i): $(max_update)/n")
    end

    println("\nSTARTING HORIZONTAL ITERATION ---")
    mesh.topleft.h, max_update =
        move_horizontally(mesh, imageBW; f = 1.0, picture_width = Caustics_Side)
    println(" max update = $(max_update)")



    println("$(mesh.topleft.vr[r_check, c_check]) , $(mesh.topleft.vc[r_check, c_check])")


    # solidMesh = create_solid(mesh)
    # save_stl!(
    #     solidMesh,
    #     "./examples/original_image.obj",
    #     scale = Float64(1 / Grid_Definition * Artifact_Size),
    #     scaleh = Float64(1 / Grid_Definition * Artifact_Size),
    # )

    return mesh, imageBW
end



"""
$(SIGNATURES)
"""
function move_vertically!(mesh, image, suffix)
    # Remember mesh is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.

    # The idea is that:
    # - any pixel on the caustic projection receives light from a given 'rectangle' on the lens.
    # - That rectangle is made of 2 triangles.
    area_distorted_pixels = get_area_pixels(mesh)
    @assert max_sum_abs(area_distorted_pixels) "ALERT: vert. Δ area_distorted_pixels is too large - max/min = $(maximum(area_distorted_pixels)) / $(minimum(area_distorted_pixels))"

    intensity_error = Float64.(area_distorted_pixels - image)
    @assert max_sum_abs(intensity_error) "ALERT: vert. Δ intensity_error is too large - max/min = $(maximum(intensity_error)) / $(minimum(intensity_error))"

    # The divergence matrix moves all the topleft corners. illumination_error is therefore not large enough.
    height, width = size(mesh)
    divergence_intensity = zeros(Float64, height + 1, width + 1)
    divergence_intensity[1:height, 1:width] .= intensity_error[1:height, 1:width]
    fill_borders!(divergence_intensity, 0.0)

    @assert max_sum_abs(divergence_intensity) "ALERT: vert. Δ divergence_intensity is too large - max/min = $(maximum(divergence_intensity)) / $(minimum(divergence_intensity))"


    # Save the loss image as a png
    plot_loss!(divergence_intensity, suffix, image)

    # ϕ is the _heightmap_ therefore size of topleft
    ϕ = mesh.topleft.h

    max_update = 100.0
    interaction_count = 0
    while true
        interaction_count += 1
        interaction_count % 500 == 0 && println("Converging for intensity: $(max_update)")

        old_max_update = max_update
        ϕ, max_update = relax_vertically(ϕ, divergence_intensity)
        if abs(max_update) < 1e-6 ||
           abs((max_update - old_max_update) / old_max_update) < 0.01
            println(
                "Convergence stopped at step $(interaction_count) with ",
                "max_update = $(max_update) -- ",
                "max/min area_distorted_pixels = $(maximum(area_distorted_pixels)) / $(minimum(area_distorted_pixels)) -- ",
                "max divergence = $(maximum(abs.(divergence_intensity)))",
            )
            break
        end
    end

    mesh.topleft.h[1:height+1, 1:width+1] .= ϕ[1:height+1, 1:width+1]

    save_stl!(
        matrix_to_mesh(ϕ * 0.02),
        "./examples/phi_$(suffix).obj",
        reverse = false,
        flipxy = true,
    )

    # plot_as_quiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reverser=false, reversec=false)
    plot_as_quiver(
        mesh,
        stride = 30,
        scale = 1.0,
        max_length = 200,
        flipxy = true,
        reverser = false,
        reversec = false,
    )

    # plot_as_quiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reversex=false, reversey=false)
    # saveObj(matrix_to_mesh(D * 10), "D_$(suffix).obj")

    # Now we need to march the x,y locations in our mesh according to this gradient!
    # Tasks are a control flow feature that allows computations to be suspended and resumed in a flexible manner.
    # We mention them here only for completeness; for a full discussion see Asynchronous Programming.
    march_mesh!(mesh, ϕ)
    save_stl!(mesh, "./examples/mesh_$(suffix).obj", flipxy = true)

    return max_update
end




"""
$(SIGNATURES)
"""
function move_horizontally(
    mesh::FaceMesh,
    image;
    f = Focal_Length,
    picture_width = Caustics_Side,
)

    # height indexes rows, width indexes columns
    height, width = size(mesh)

    # Coordinates difference
    # Coordinates on the caustics = simple rectangular values in pixels.
    # Coordinates on the lens face comes from topleft field.
    d_row = mesh.topleft.rows_numbers - mesh.topleft.r
    d_col = mesh.topleft.cols_numbers - mesh.topleft.c

    # true_H is the real distance beteen the surface of the lens and the projection
    # H is the initial distance from the surface of the carved face to the projection screen.
    # H is constant and does not account for the change in heights due to carving.
    # The focal length is in meters. H is in Pixels.
    # Note: Higher location means closer to the caustics => negative sign
    H = f / Meters_Per_Pixel

    true_H = similar(mesh.topleft.h)
    @. true_H[:] = -mesh.topleft.h[:] + H

    # Normals.
    N_row = zeros(Float64, height + 1, width + 1)
    N_col = zeros(Float64, height + 1, width + 1)
    @. N_row[:] = tan(atan(d_row[:] / true_H[:]) / (n₁ - n₂))
    @. N_col[:] = tan(atan(d_col[:] / true_H[:]) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    div_row = zeros(Float64, width, height)
    div_col = zeros(Float64, width, height)
    @. div_row[1:end, 1:end] = N_row[2:end, 1:end-1] - N_row[1:end-1, 1:end-1]
    @. div_col[1:end, 1:end] = N_col[1:end-1, 2:end] - N_col[1:end-1, 1:end-1]

    divergence_direction = zeros(Float64, width + 1, height + 1)
    divergence_direction[1:height, 1:width] = div_row + div_col
    fill_borders!(divergence_direction, 0.0)

    @assert max_sum_abs(divergence_direction) "ALERT: horiz. Δ divergence_direction is too large - max/min = $(maximum(divergence_direction)) / $(minimum(divergence_direction))"


    println("Have all the divergences")

    ϕ = zeros(Float64, height + 1, width + 1)
    max_update = 1_000.0
    iteration_count = 0
    while true
        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for divergence: $(max_update)")

        old_max_update = max_update
        ϕ, max_update = relax_vertically(ϕ, -divergence_direction)
        @assert max_sum_abs(ϕ) "ALERT: horiz. Δ ϕ is too large - max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))"


        if abs(max_update) < 1e-6 ||
           abs((max_update - old_max_update) / old_max_update) < 0.01
            println(
                "Convergence stopped improving at step $(iteration_count) ",
                "with max_update of $(max_update)",
            )
            break
        end
    end

    # saveObj(matrix_to_mesh(h / 10), "./examples/heightmap.obj")
    return ϕ, max_update
end


"""
$(SIGNATURES)

This function implements successive over relaxation for a matrix and its associated error matrix
There is a hardcoded assumption of Neumann boundary conditions--that the derivative across the
boundary must be zero in all cases. See:
https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

"""
function relax_vertically(ϕ::Matrix{Float64}, divergence::Matrix{Float64})

    @assert max_sum_abs(ϕ) "ALERT: relaxation ϕ is too large - max/min = $(maximum(ϕ)) / $(minimum(ϕ))"
    @assert max_sum_abs(divergence) "ALERT: relaxation divergence is too large - max/min = $(maximum(divergence)) / $(minimum(divergence))"

    @assert any(isnan.(ϕ)) == false "Calling relaxation with a ϕ that contains NaN!!"
    @assert any(isnan.(divergence)) == false "Calling relaxation with a divergence that contains NaN!!"


    # ϕ is of the same size of topleft
    height, width = size(ϕ)
    n_pixels = height * width
    height_average = sum(ϕ) / n_pixels

    # Laplacian

    # Embed matrix within a larger matrix for better vectorization and avoid duplicated code
    # Size of ϕ/topleft is h+1, w+1. container adds 1 row/col before and 1 row/col after
    container_height = 1 + height + 1
    container_width = 1 + width + 1
    container = zeros(Float64, container_height, container_width)
    container[2:container_height-1, 2:container_width-1] .= ϕ
    @assert max_sum_abs(container) "ALERT: relaxation container is too large - max/min = $(maximum(container)) / $(minimum(container))"

    # Those values are all of ϕ's size but shifted in 4 directions
    height_above = container[2-1:container_height-1-1, 2+0:container_width-1]
    height_below = container[2+1:container_height-1+1, 2+0:container_width-1]
    height_left = container[2+0:container_height-1, 2-1:container_width-1-1]
    height_right = container[2+0:container_height-1, 2+1:container_width-1+1]

    # Target position. The target is the current height map smoothed by averaging to which the
    # divergence is added.
    # δ_map is the same size as divergence.
    δ_map = zeros(Float64, size(divergence))

    # δ_map is the Laplacian of the height
    δ_map = (height_above + height_below + height_left + height_right) / 4.0 - ϕ
    @assert max_sum_abs(δ_map) "ALERT: relaxation δ_map as (average - ϕ) is too large - max/min = $(maximum(δ_map)) / $(minimum(δ_map))"

    # This has to converge towards the divergence. Difference to calculate speed of the descent.
    δ_map -= divergence
    @assert max_sum_abs(δ_map) "ALERT: relaxation δ_map as including divergence is too large - max/min = $(maximum(δ_map)) / $(minimum(δ_map))"

    # height_div, width_div = size(divergence)
    # @. target_map[height_div, width_div] += divergence / 4.0
    fill_borders!(δ_map, 0.0)

    # Let the heightmap converge towards the target at a slow rate.
    # new_ϕ = ϕ + ω .* δ_map
    new_ϕ = zeros(Float64, size(divergence))

    # new_ϕ = ϕ + 1.94 * δ_map
    new_ϕ = ϕ + 0.2 * δ_map
    @assert any(isnan.(new_ϕ)) == false "new_ϕ contains NaN!!"

    max_update = maximum(abs.(new_ϕ))
    @assert !isnan(max_update) """ MAX UPDATE IS NOT A NaN!!!
                                   Max height = $( maximum(ϕ) )  --  Min height = $( minimum(ϕ) )  --  Height average = $(height_average)
                                   Max divergence = $( maximum(divergence) )  --  Min divergence = $( minimum(divergence) )
                                   """

    return new_ϕ, max_update
end


"""
$(SIGNATURES)

A Mesh is a collection of triangles. The brightness flowing through a given triangle is just proportional to its
area in the x, y plane. h is ignored.

The function returns a matrix with the quantity of light coming from each 'rectangle'  around a pixel. That 'rectangle'
has been shifted and flexed around.
"""
function get_area_pixels(mesh::FaceMesh)
    height, width = size(mesh)

    pixel_areas = zeros(Float64, size(mesh))

    for ci ∈ CartesianIndices(pixel_areas)
        top_tri_area = area(triangle3D(mesh, ci, :top)...)
        bot_tri_area = area(triangle3D(mesh, ci, :bottom)...)

        isnan(top_tri_area) &&
            println("NaN area at $(ci) for top triangle $(top_triangle3D(mesh, ci)).")
        isnan(bot_tri_area) &&
            println("NaN area at $(ci) for bottom triangle $(bot_triangle3D(mesh, ci)).")

        pixel_areas[ci] = top_tri_area + bot_tri_area
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

    normalised_D_max = D ./ maximum(D)
    normalised_D_min = D ./ minimum(D)

    blue = zeros(size(D))
    blue[D.>0] = normalised_D_max[D.>0]
    red = zeros(size(D))
    red[D.<0] = -normalised_D_min[D.<0]
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
    Br = p2.r - p1.r
    Bc = p2.c - p1.c
    Cr = p3.r - p1.r
    Cc = p3.c - p1.c

    t_vBr = p2.vr - p1.vr
    t_vBc = p2.vc - p1.vc
    t_vCr = p3.vr - p1.vr
    t_vCc = p3.vc - p1.vc

    # After this, given that Ar = Ac = t_vAr = t_vAc = 0, the area is nil iff
    # (Br + t_vBr) (Cc + t_vCc ) - (Cr + t_vCr) (Bc + t_vBc) = 0.
    # After expansion and reshuffling to have a quadratic equation where t
    # is the variable, the coefficients of that equation are:
    a = t_vCc * t_vBr - t_vBc * t_vCr
    b = -Bc * t_vCr - Cr * t_vBc + Br * t_vCc + Cc * t_vBr
    c = Br * Cc - Cr * Bc

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
            # C will never end up on the line AB.
            # Very unlikely with Float64.
            # Negative numbers are filtered out when calculating the minimum jiggle ratio
            return -1.0, -1.0
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
    # println("Min/Max of values of ϕ: $(minimum(ϕ)) / $(maximum(ϕ))")
    # println("Sum of values of ∇ϕᵤ:   $( sum(∇ϕᵤ) )")
    # println("Sum of values of ∇ϕᵥ:   $( sum(∇ϕᵥ) )")

    height, width = size(ϕ)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.
    # Warning: The indices are necessary because topleft and ∇ϕ are of different sizes
    # CHECK SIGNS!
    mesh.topleft.vr[1:height, 1:width] .= -∇ϕᵤ[1:height, 1:width]
    mesh.topleft.vc[1:height, 1:width] .= -∇ϕᵥ[1:height, 1:width]

    # Just in case...
    reset_border_values!(mesh.topleft)

    # Basically infinity time to shrink a triangle to nothing.
    min_positive_t = Inf
    list_min_pos_t = zeros(Float64, (height - 2) * (width - 2))

    t1 = zeros(Float64, height, width)
    t2 = zeros(Float64, height, width)

    mesh_r = copy(mesh.topleft.r)
    mesh_c = copy(mesh.topleft.c)
    mesh_h = copy(mesh.topleft.h)
    for row ∈ 1:height-1, col ∈ 1:width-1
        # # Jiggle each point in a random order.
        # for (row, col) ∈ shuffle([(r, c) for r ∈ 2:height-1, c ∈ 2:width-1])

        # Get the time, at that velocity, for the area of the triangle to be nil.
        # We are only interested by positive values to only move in the direction of the gradient
        for triangle ∈
            [triangle3D(mesh, row, col, :top), triangle3D(mesh, row, col, :bottom)]
            t1, t2 = find_maximum_t(triangle)
            @assert typeof(t1) == Float64 && typeof(t2) == Float64 "Maximum times are not numerical at $(row), $(col)"
        end

        # list_min_pos_t[(row-1)+(col-2)*(width-2)]
        # println("March mesh δ: $(δ)")
    end

    t1 = max.(t1, 0.0)
    t2 = max.(t2, 0.0)

    min_positive_t = min(minimum(t1), minimum(t2))

    δ = min_positive_t / 2.0
    mesh.topleft.r -= δ .* ∇ϕᵤ
    mesh.topleft.c -= δ .* ∇ϕᵥ

    # reset_border_values!(mesh.topleft)

    # Reset the border at the fixed values fixed coordinates.
    reset_border_values!(mesh.topleft)

    println(
        "Average mesh changes on row = $(sum(abs.(mesh.topleft.r - mesh_r)) / (height * width))",
    )

    println(
        "Average mesh changes on col = $(sum(abs.(mesh.topleft.c - mesh_c)) / (height * width))",
    )

    println(
        "Average mesh changes on height = $(sum(abs.(mesh.topleft.h - mesh_h)) / (height * width))",
    )

    # Modify the mesh triangles but ensuring that we are far from destroying any of them with a nil area.
    # Use the maximum change possible (negative or positive)
    # δ = min_positive_t / 2.0
    # println(
    #     "\tMarch mesh with average factor $(sum(list_min_pos_t) / length(list_min_pos_t))\n",
    # )

    # println("Overall maximum variations:")
    # println("\tOverall min_positive_t: $(min_positive_t)")
    # println("\tOverall max_negative_t: $(max_negative_t)")
    # println("\tUsing δ: $(δ)")

    # saveObj(mesh, "gateau.obj")
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
function matrix_to_mesh(ϕ::Matrix{Float64})
    height, width = size(ϕ)

    mesh = FaceMesh(height, width)

    # 1 more topleft than pixels! Therefore compiler needs to specify exact indices.
    mesh.topleft.h[1:height, 1:width] .= ϕ[1:height, 1:width]

    # The borders' height is forced at 0.
    reset_border_values!(mesh.topleft)

    return mesh
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
    bottom_mesh.topleft.r[:] .= rows[:]
    bottom_mesh.topleft.c[:] .= cols[:]
    bottom_mesh.topleft.h[:] .= -Bottom_Offset
    bottom_mesh.topleft.vr[:] .= 0.0
    bottom_mesh.topleft.vr[:] .= 0.0


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
