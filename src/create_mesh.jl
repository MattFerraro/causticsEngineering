"""
$(SIGNATURES)
"""
function engineer_caustics(source_image)
    imageBW = Float64.(Gray.(source_image))

    height, width = size(imageBW)
    println("Image size: $((height, width))")

    # mesh is the same size as the image with an extra row/column to have coordinates to
    # cover each image corner with a triangle.
    mesh = FaceMesh(height, width)
    # println("Mesh creation: ")
    # r_check = 1
    # c_check = 1
    # println(
    #     "\tCoordinates: $(mesh.corners.r[r_check, c_check]) , $(mesh.corners.c[r_check, c_check]) , $(mesh.corners.ϕ[r_check, c_check])",
    # )
    # println(
    #     "\tVelocities: $(mesh.corners.vr[r_check, c_check]) , $(mesh.corners.vc[r_check, c_check])",
    # )

    # The energy going through the lens is equal to the amount of energy on the caustics
    total_energy_lens = height * width * 1 # 1 unit of energy per corner
    total_energy_caustics = sum(imageBW)
    correction_ratio = (width * height) / sum(imageBW)

    # imageBW is `grid_definition x grid_definition` and is normalised to the same (sort of) _energy_ as the
    # original image.
    imageBW = imageBW .* correction_ratio

    for i ∈ 1:5
        println("\nSTARTING VERTICAL ITERATION $(i) ---")

        max_update = solve_velocity_potential!(mesh, imageBW, "it$(i)")
        print("Vertical move max update = $(max_update)")

        # println(
        #     "\tCoordinates: $(mesh.corners.r[r_check, c_check]) , $(mesh.corners.c[r_check, c_check]) , $(mesh.corners.ϕ[r_check, c_check])",
        # )
        # println(
        #     "\tVelocities: $(mesh.corners.vr[r_check, c_check]) , $(mesh.corners.vc[r_check, c_check])",
        # )

        println("---------- ITERATION $(i): $(max_update)/n")
    end

    println("\nSTARTING HORIZONTAL ITERATION ---")
    mesh.corners.ϕ, max_update =
        move_horizontally(mesh, imageBW; f = 1.0, picture_width = Caustics_Side)
    println(" max update = $(max_update)")



    println("$(mesh.corners.vr[r_check, c_check]) , $(mesh.corners.vc[r_check, c_check])")


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
function solve_velocity_potential!(mesh, image, suffix)
    # Remember mesh is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.

    # The idea is that:
    # - any corner on the caustic projection receives light from a given 'rectangle' on the lens.
    # - That rectangle is made of 2 triangles.
    area_distorted_corners = get_area_corners(mesh)
    @assert average_absolute(area_distorted_corners) """

                                                    solve_velocity_potential: area_distorted_corners values seem too large
                                                    $(maximum(area_distorted_corners)) / $(minimum(area_distorted_corners))

                                                        """


    error_luminosity = Float64.(area_distorted_corners - image)
    @assert average_absolute(error_luminosity) """

                                                   solve_velocity_potential: error_luminosity values seem too large
                                                       max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))
                                                       max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

                                                       """

    # Save the loss image as a png
    plot_loss!(error_luminosity, suffix, image)

    # ϕ is the _velocity potential_
    # ϕ is defined per corner
    ϕ = zeros(Float64, size(mesh.corners.ϕ))

    # For the purpose of Poisson, we need to divergence to increase towards the inside: negative divergence for high luminosity
    error_luminosity = -error_luminosity

    max_update = Inf
    iteration_count = 0
    for _ ∈ 1:5
        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for intensity: $(max_update)")

        old_max_update = max_update
        println(
            "solve_velocity_potential: iteration count = $(iteration_count) -- ",
            "max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ)) -- ",
            "max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ)) -- ",
            "max_update = $(max_update)",
        )

        # Solve the Poisson equation where the divergence of the gradient of ϕ is equal to the luminosity loss.
        ϕ, max_update = gradient_descent(ϕ, error_luminosity)
        @assert average_absolute(ϕ) """

                                    solve_velocity_potential: ϕ values seem too large
                                        iteration count = $(iteration_count)
                                        max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))
                                        max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))

                                        max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))
                                        max_update = $(max_update)

                                        """



        if abs(max_update) < 1e-6 ||
           abs((max_update - old_max_update) / old_max_update) < 0.01
            println(
                "Convergence stopped at step $(iteraction_count) with ",
                "max_update = $(max_update) -- ",
                "max/min area_distorted_corners = $(maximum(area_distorted_corners)) / $(minimum(area_distorted_corners)) -- ",
                "max error_luminosity = $(maximum(abs.(error_luminosity)))",
            )
            break
        end
    end

    mesh.corners.ϕ .= ϕ

    # Now we need to march the x,y locations in our mesh according to this gradient!
    # Tasks are a control flow feature that allows computations to be suspended and resumed in a flexible manner.
    # We mention them here only for completeness; for a full discussion see Asynchronous Programming.
    march_mesh!(mesh)

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
    # Coordinates on the caustics = simple rectangular values in corners.
    # Coordinates on the lens face comes from corner field.
    d_row = mesh.corners.rows_numbers - mesh.corners.r
    d_col = mesh.corners.cols_numbers - mesh.corners.c

    # true_H is the real distance beteen the surface of the lens and the projection
    # H is the initial distance from the surface of the carved face to the projection screen.
    # H is constant and does not account for the change in heights due to carving.
    # The focal length is in meters. H is in corners.
    # Note: Higher location means closer to the caustics => negative sign
    H = f / Meters_Per_corner

    true_H = similar(mesh.corners.ϕ)
    @. true_H[:] = -mesh.corners.ϕ[:] + H

    # Normals.
    N_row = zeros(Float64, size(true_H))
    N_col = zeros(Float64, size(true_H))
    @. N_row[:] = tan(atan(d_row[:] / true_H[:]) / (n₁ - n₂))
    @. N_col[:] = tan(atan(d_col[:] / true_H[:]) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    div_row = zeros(Float64, width, height)
    div_col = zeros(Float64, width, height)
    @. div_row[1:end, 1:end] = N_row[2:end, 1:end-1] - N_row[1:end-1, 1:end-1]
    @. div_col[1:end, 1:end] = N_col[1:end-1, 2:end] - N_col[1:end-1, 1:end-1]

    divergence_direction = zeros(Float64, size(div_row))
    divergence_direction[1:height, 1:width] = div_row + div_col
    @assert average_absolute(divergence_direction) "ALERT: horiz. Δ divergence_direction is too large - max/min = $(maximum(divergence_direction)) / $(minimum(divergence_direction))"

    println("Have all the divergences")

    corner_heights = zeros(Float64, size(mesh.corners.ϕ))
    max_update = 1_000.0
    iteration_count = 0
    while true
        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for divergence: $(max_update)")

        old_max_update = max_update
        corner_heights, max_update = gradient_descent(corner_heights, -divergence_direction)
        @assert average_absolute(corner_heights) "ALERT: horiz. Δ corner_heights is too large - max/min corner_heights = $(maximum(corner_heights)) / $(minimum(corner_heights))"


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
    return corner_heights, max_update
end


"""
$(SIGNATURES)

This function implements successive over relaxation for a matrix and its associated error matrix
There is a hardcoded assumption of Neumann boundary conditions--that the derivative across the
boundary must be zero in all cases. See:
https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

- ϕ is the potential to solve subject to the Poisson equation and is the same size as the corners (posts)
- target is ∇²ϕ and is of the size of the pixels (fences)

"""
function gradient_descent(ϕ::Matrix{Float64}, ∇ϕ::Matrix{Float64})

    # Ensure that border conditions are as they should
    fill_borders!(∇ϕ, 0.0)
    fill_borders!(ϕ, 0.0)

    @assert average_absolute(ϕ) """

                                Relaxation: ϕ values seem too large
                                    max/min ∇ϕ = $(maximum(∇ϕ)) / $(minimum(∇ϕ))
                                    max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

                                    """

    @assert average_absolute(∇ϕ) """

                                         Relaxation: ∇²ϕ values seem too large
                                             max/min ∇ϕ = $(maximum(∇ϕ)) / $(minimum(∇ϕ))
                                             max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

                                             """


    @assert any(isnan.(ϕ)) == false "Relaxation: calling  with a ϕ that contains NaN!!"
    @assert any(isnan.(∇ϕ)) == false "Relaxation: calling with a ∇ϕ that contains NaN!!"


    # ϕ is of the same size as corners = number of corners + 1
    height_corners, width_corners = size(ϕ)
    height_pixels, width_pixels = size(∇ϕ)

    # Laplacian

    # Embed matrix within a larger matrix for better vectorization and avoid duplicated code
    # Size of ϕ/corner is h, w. Padded matrix adds 1 row/col after
    # ϕ is inserted in padded matrix within  1:height_corners+1 x 1:width_corners+1.
    # The rest of the padded matrix (the borders) are set at 0.0.
    padded_ϕ = zeros(Float64, height_corners + 1, width_corners + 1)
    padded_ϕ[1:height_corners, 1:width_corners] .= ϕ

    # Those values are all of ϕ's size and represent the _flow_ in each direction
    flow_up = zeros(Float64, size(∇ϕ))
    flow_down = zeros(Float64, size(∇ϕ))
    flow_left = zeros(Float64, size(∇ϕ))
    flow_right = zeros(Float64, size(∇ϕ))

    flow_up[1:end, 1:end] .= padded_ϕ[1:end-2, 2:end-1] - padded_ϕ[2:end-1, 2:end-1]
    flow_down[1:end, 1:end] .= padded_ϕ[3:end, 2:end-1] - padded_ϕ[2:end-1, 2:end-1]
    flow_left[1:end, 1:end] .= padded_ϕ[2:end-1, 1:end-2] - padded_ϕ[2:end-1, 2:end-1]
    flow_right[1:end, 1:end] .= padded_ϕ[2:end-1, 3:end] - padded_ϕ[2:end-1, 2:end-1]

    fill_borders!(flow_up, 0.0)
    fill_borders!(flow_down, 0.0)
    fill_borders!(flow_left, 0.0)
    fill_borders!(flow_right, 0.0)

    # Target position. The target is the current height map smoothed by averaging to which the
    # flow is added.
    # This has to converge towards the ∇ϕ. Difference to calculate speed of the descent.
    δ = (flow_up + flow_down + flow_left + flow_right) + ∇ϕ

    @assert average_absolute(δ) """

                                               Relaxation: ∇ϕ_error as (current - target) seem too large
                                                   max/min target ∇ϕ = $(maximum(∇ϕ)) / $(minimum(∇ϕ))
                                                   max/min ∇ϕ_error = $(maximum(δ)) / $(minimum(δ))

                                                   max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))
                                                   max/min padded_ϕ = $(maximum(padded_ϕ)) / $(minimum(padded_ϕ))
                                                   max/min flow_up = $(maximum(flow_up)) / $(minimum(flow_up))
                                                   max/min flow_down = $(maximum(flow_down)) / $(minimum(flow_down))
                                                   max/min flow_left = $(maximum(flow_left)) / $(minimum(flow_left))
                                                   max/min flow_right = $(maximum(flow_right)) / $(minimum(flow_right))

                                                   """

    # @. target_map[height_div, width_div] += ∇ϕ / 4.0
    # fill_borders!(∇ϕ_error, 0.0)

    # Let the heightmap converge towards the target at a slow rate = gradient descent.
    # @. ϕ += 0.2 * ∇ϕ_error
    @. ϕ[1:height_pixels, 1:width_pixels] +=
        1.94 / 40.0 * δ[1:height_pixels, 1:width_pixels]
    @assert !any(isnan.(ϕ)) "Relaxation: Updated ϕ contains NaN!!"

    max_update = maximum(abs.(1.94 / 40.0 * δ))

    return ϕ, max_update
end


"""
$(SIGNATURES)

A Mesh is a collection of triangles. The brightness flowing through a given triangle is just proportional to its
area in the x, y plane. h is ignored.

The function returns a matrix with the quantity of light coming from each 'rectangle'  around a corner. That 'rectangle'
has been shifted and flexed around.
"""
function get_area_corners(mesh::FaceMesh)
    height, width = size(mesh)

    top_tri_area =
        [area(triangle3D(mesh, row, col, :top)...) for row ∈ 1:height, col ∈ 1:width]
    @assert !any(isnan.(top_tri_area)) "get_area_corners: NaN area in top triangles."

    bot_tri_area =
        [area(triangle3D(mesh, row, col, :bottom)...) for row ∈ 1:height, col ∈ 1:width]
    @assert !any(isnan.(top_tri_area)) "get_area_corners: NaN area in bottom triangles."

    return top_tri_area + bot_tri_area
end


"""
$(SIGNATURES)
"""
function quantify_loss(D, suffix, img)
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

    rgbImg = RGB.(red, green, blue)'
    save("./examples/loss_$(suffix).png", map(clamp01nan, rgbImg))

    println("Saving output image:")
    println(typeof(img))
    E = Gray.(D)
    println(typeof(E))
    outputImg = img - E
    save("./examples/actual_$(suffix).png", outputImg)
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
    if a == 0 && b != 0
        return smallest_positive(-c / b, c / b)
    else
        discriminant = b^2 - 4a * c

        # If there is a solution
        if discriminant >= 0
            d = sqrt(discriminant)
            return smallest_positive((-b - d) / 2a, (-b + d) / 2a)
        end
    end
    # There can be no solution if, after translation, B abd C move in parallel direction.
    # C will never end up on the line AB.
    # Very unlikely with Float64.
    # Negative numbers are filtered out when calculating the minimum jiggle ratio
    return 1e-10
end


find_maximum_t(p::Tuple{Vertex3D,Vertex3D,Vertex3D}) = find_maximum_t(p[1], p[2], p[3])


"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh!(mesh::FaceMesh)

    ϕ = mesh.corners.ϕ
    height, width = size(ϕ)

    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)
    # println("Min/Max of values of ϕ: $(minimum(ϕ)) / $(maximum(ϕ))")
    # println("Sum of values of ∇ϕᵤ:   $( sum(∇ϕᵤ) )")
    # println("Sum of values of ∇ϕᵥ:   $( sum(∇ϕᵥ) )")


    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.
    # Warning: The indices are necessary because corner and ∇ϕ are of different sizes
    # CHECK SIGNS!
    mesh.corners.vr .= -∇ϕᵤ
    mesh.corners.vc .= -∇ϕᵥ

    # Just in case...
    reset_border_values!(mesh.corners)
    reset_border_values!(mesh.corners)

    # Basically infinity time to shrink a triangle to nothing.
    min_positive_t = Inf
    list_min_pos_t = zeros(Float64, (height - 2) * (width - 2))

    t1 = zeros(Float64, size(mesh))
    t2 = zeros(Float64, size(mesh))

    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)
    mesh_h = copy(mesh.corners.ϕ)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the gradient

    list_triangles = vcat(
        [triangle3D(mesh, row, col, :top) for row ∈ 1:height-1, col ∈ 1:width-1],
        [triangle3D(mesh, row, col, :bottom) for row ∈ 1:height-1, col ∈ 1:width-1],
    )
    list_maximum_t = find_maximum_t.(list_triangles)
    min_positive_t = minimum(list_maximum_t)

    # @assert all(typeof.(t1) .== Float64) && all(typeof.(t2) .== Float64) "March mesh: Maximum times are not numerical at $(row), $(col)"

    δ = min_positive_t / 2.0
    mesh.corners.r -= δ .* ∇ϕᵤ
    mesh.corners.c -= δ .* ∇ϕᵥ

    # reset_border_values!(mesh.corner)

    # Reset the border at the fixed values fixed coordinates.
    reset_border_values!(mesh.corners)

    println(
        """

        March mash:
            δ = $(δ)
            first minimum ts: $(list_maximum_t[1:5])
            Average mesh changes on row = $(sum(abs.(mesh.corners.r - mesh_r)) / (height * width))
            Average mesh changes on col = $(sum(abs.(mesh.corners.c - mesh_c)) / (height * width))
            Average mesh changes on height = $(sum(abs.(mesh.corners.ϕ - mesh_h)) / (height * width))

        """,
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
    ∇ϕᵤ = zeros(Float64, size(ϕ))   # divergence on the right edge will be filled with zeros
    ∇ϕᵥ = zeros(Float64, size(ϕ))   # divergence on bottom edge will be filled with zeros

    @. ∇ϕᵤ[begin:end-1, :] = ϕ[begin+1:end, :] - ϕ[begin:end-1, :]
    @. ∇ϕᵥ[:, begin:end-1] = ϕ[:, begin+1:end] - ϕ[:, begin:end-1]

    fill_borders!(∇ϕᵤ, 0.0)
    fill_borders!(∇ϕᵥ, 0.0)

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

    # 1 more corner than corners! Therefore compiler needs to specify exact indices.
    mesh.corners.ϕ[1:height, 1:width] .= ϕ[1:height, 1:width]

    # The borders' height is forced at 0.
    reset_border_values!(mesh.corners)

    return mesh
end


"""
$(SIGNATURES)

TO REFACTOR.
"""
function create_solid(
    mesh::FaceMesh;
    bottom_offset = Bottom_Offset / Meters_Per_corner,
    top_offset = Top_Offset / Meters_Per_corner,
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
    bottom_mesh.corner.r[:] .= rows[:]
    bottom_mesh.corner.c[:] .= cols[:]
    bottom_mesh.corner.ϕ[:] .= -Bottom_Offset
    bottom_mesh.corner.vr[:] .= 0.0
    bottom_mesh.corner.vr[:] .= 0.0


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
