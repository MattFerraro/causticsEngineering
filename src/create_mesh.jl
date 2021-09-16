
"""
$(SIGNATURES)

"""
function engineer_caustics(source_image; clamp_correction = true)
    imageBW = Float64.(Gray.(source_image))


    # Set global size parameters
    height, width = size(imageBW)
    global N_Pixel_Height = height
    global N_Pixel_Width = width
    global Meters_Per_Pixel = Caustics_Long_Side / N_Pixel_Height
    println("Image size: $(N_Pixel_Height) x $(N_Pixel_Width)")


    # mesh is the same size as the image with an extra row/column to have coordinates to
    # cover each image corner with a triangle.
    mesh = FaceMesh(N_Pixel_Height, N_Pixel_Width)

    # The total energy going through the lens is equal to the amount of energy on the caustics
    # imageBW is normalised to the same (sort of) _energy_ as the original image.
    imageBW = imageBW / average(imageBW)

    # start_max_update has no particular meaning.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update

    error_luminosity = 1e6 * N_Pixel_Height * N_Pixel_Width
    old_error_luminosity = 2 * error_luminosity
    iteration_count = 0

    while (1e-6 < abs(max_update) < 10 * start_max_update &&
           # 1e-4 < abs((max_update - old_max_update) / old_max_update) &&
           iteration_count < 1_024)
        iteration_count += 1

        println(
            """
            engineer_caustics:
                STARTING VERTICAL ITERATION $(iteration_count) --------------------------------------------

                """,
        )

        old_max_update = max_update
        old_error_luminosity = error_luminosity
        max_update, error_luminosity = solve_velocity_potential!(
            mesh,
            imageBW,
            "loss_it$(iteration_count)";
            clamp_correction = clamp_correction,
        )
        print(
            """
              engineer_caustics:
                  Vertical move max update = $(max_update)   --  Mean error luminosity = $(error_luminosity)
                  """,
        )
    end

    println("\nSTARTING HORIZONTAL ITERATION AFTER  $(iteration_count) iterations ---\n")

    ϕ, max_update =
        move_horizontally(mesh, imageBW; f = Focal_Length, clamp_correction = true)
    mesh.corners.ϕ .= ϕ

    println("""
    engineer_caustics:
        max/min ϕ = $(maximum(abs.(ϕ))) / $(minimum(abs.(ϕ)))
        size ϕ = $(size(ϕ))
        average ϕ = $(average_absolute(ϕ))
        Horizontal max update = $(max_update)
        """)

    # Move the around a nil average.
    mean_ϕ = average(mesh.corners.ϕ)
    mesh.corners.ϕ = mesh.corners.ϕ .- mean_ϕ

    return mesh, imageBW
end



"""
$(SIGNATURES)

ϕ is the _velocity potential_. The velocity represents the direction towards which the mesh vertices should move.
Zones of high luminosity should be covered with more triangles and should attract more vertices.
"""
function solve_velocity_potential!(mesh, image, prefix; clamp_correction = true)
    # Keep that to track how much progress is made
    mesh_ϕ = copy(mesh.corners.ϕ)

    # Get the area of each individual pixel as stretch/shrunk on the lens. Area = energy.
    # _FENCES_SIZED_
    lens_pixels_area = get_lens_pixels_area(mesh)

    # Positive error => the triangle needs to shrink (less light)
    error_luminosity = Float64.(lens_pixels_area - image)
    println(
        """
            solve_velocity_potential:
                max/min error_luminosity: $(maximum(error_luminosity)) / $(minimum(error_luminosity))
                """,
    )
    @assert test_average_absolute(error_luminosity) """

            solve_velocity_potential: error_luminosity values seem too large
                max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))

            """
    save_plot_scalar_field!(error_luminosity, prefix * "_loss", image)


    # start_max_update has no particular meeting.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update
    iteration_count = 0
    new_divergence = 0.0

    while (
        1e-5 < abs(max_update) < 10 * start_max_update &&
        1e-2 < abs((max_update - old_max_update) / old_max_update) &&
        new_divergence < 100.0
    )

        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for intensity: $(max_update)")

        old_max_update = max_update
        println("""
            solve_velocity_potential: before propagation
                max/min ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))
                iteration count = $(iteration_count)
                max_update = $(max_update)
            """)

        # Solve the Poisson equation where the divergence of the gradient of ϕ is equal to the luminosity error.
        # High error => Too much luminosity.
        mean_ϕ = sum(mesh.corners.ϕ) / length(mesh.corners.ϕ)
        mesh.corners.ϕ .-= mean_ϕ

        ϕ, max_update = propagate_poisson(
            mesh.corners.ϕ,
            error_luminosity;
            clamp_correction = clamp_correction,
        )

        mean_ϕ = sum(ϕ) / length(ϕ)
        mesh.corners.ϕ = ϕ .- mean_ϕ


        println(
            """
        solve_velocity_potential: after propagation
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))
            iteration count = $(iteration_count)
            max_update = $(max_update)
        """,
        )

        @assert test_average_absolute(mesh.corners.ϕ) """
                solve_velocity_potential: ϕ values seem too large
                    iteration count = $(iteration_count)
                    max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))

                    max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))
                    max_update = $(max_update)

                    """
        new_divergence = maximum(
            abs.(
                mesh.corners.ϕ[begin+1:end, begin:end-1] .-
                mesh.corners.ϕ[begin:end-1, begin:end-1] .+
                mesh.corners.ϕ[begin:end-1, begin+1:end] .-
                mesh.corners.ϕ[begin:end-1, begin:end-1],
            ),
        )

    end

    println(
        """

    Luminosity convergence stopped at step $(iteration_count) with
        max_update = $(max_update)
        max / min error_luminosity = $(maximum(abs.(error_luminosity))) / $(minimum(abs.(error_luminosity)))
        Average mesh changes on ϕ = $(average_absolute(mesh.corners.ϕ - mesh_ϕ))
        max ϕ divergence = $(maximum(abs.(mesh.corners.ϕ[begin+1:end, begin:end-1] .- mesh_ϕ[begin:end-1, begin:end-1] .+
        mesh.corners.ϕ[begin:end-1, begin+1:end] .- mesh_ϕ[begin:end-1, begin:end-1])))

        """,
    )

    # Now we need to march the mesh row,col corner locations according to this gradient.
    march_mesh!(mesh)
    plot_as_quiver(mesh, n_steps = 50, scale = 1.0, max_length = N_Pixel_Height / 20)

    return max_update, average_absolute(error_luminosity)
end


"""
$(SIGNATURES)
"""
function move_horizontally(mesh::FaceMesh, image; f = Focal_Length, clamp_correction = true)

    # height indexes rows, width indexes columns
    height, width = size(mesh)

    # Coordinates on the caustics = simple rectangular values in corners.
    # Coordinates on the lens face comes from corner field.
    # d_row/d_col are _POSTS-SIZED_
    d_row = mesh.corners.rows_numbers - mesh.corners.r
    d_col = mesh.corners.cols_numbers - mesh.corners.c

    # The focal length is in meters. H is in corners.
    H = f / Meters_Per_Pixel

    # true_H is the real distance beteen the surface of the lens and the projection
    # H is the initial distance from the surface of the carved face to the projection screen.
    # H is constant and does not account for the change in heights due to carving.
    # Note: Higher location means closer to the caustics => negative sign
    # true_H is _POSTS_SIZED_
    true_H = zeros(Float64, size(mesh.corners.ϕ))
    @. true_H = -mesh.corners.ϕ + H

    # Normals.
    # N_row/N_col are _POSTS_SIZED_
    N_row = zeros(Float64, size(true_H))
    N_col = zeros(Float64, size(true_H))
    @. N_row[1:end, 1:end] =
        tan(atan(d_row[1:end, 1:end] / true_H[1:end, 1:end]) / (n₁ - n₂))
    @. N_col[1:end, 1:end] =
        tan(atan(d_col[1:end, 1:end] / true_H[1:end, 1:end]) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    # div_row/div_col are _FENCES_SIZED_
    div_row = zeros(Float64, height, width)
    div_col = zeros(Float64, height, width)
    div_row[1:end, 1:end] = N_row[2:end, 1:end-1] .- N_row[1:end-1, 1:end-1]
    div_col[1:end, 1:end] = N_col[1:end-1, 2:end] .- N_col[1:end-1, 1:end-1]

    # divergence_direction = zeros(Float64, size(div_row))
    # divergence_direction is _FENCES_SIZED_
    divergence_direction = div_row + div_col
    @assert test_average_absolute(divergence_direction) """

            ALERT: horiz.
                Δ divergence_direction seem too large
                max/min = $(maximum(divergence_direction)) / $(minimum(divergence_direction))

                """

    ##
    ## CHECK SIGN!!!
    ##
    true_H, max_update = propagate_poisson(
        true_H,
        -divergence_direction;
        clamp_correction = clamp_correction,
    )
    @assert test_average_absolute(true_H) """

        ALERT:
            horiz. Δ corner_heights is too large
            max/min corner_heights = $(maximum(true_H)) / $(minimum(true_H))

            """

    println("Convergence horiz. Δ stopped with max_update of $(max_update)")

    return true_H, max_update
end


"""
$(SIGNATURES)

This function implements successive over relaxation for a matrix and its associated error matrix
There is a hardcoded assumption of Neumann boundary conditions -- that the derivative across the
boundary must be zero in all cases. See:
https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

- ϕ is the potential to be solved subject to the Poisson equation and is the same size as the corners (_POSTS_SIZED_)
- target is ∇ϕ and is of the size in pixels (_FENCES_SIZED_)

In order of how to think about the flow of the program:

- High ∇ϕ (intensity loss function for the velocity potential) means that the triangles are too bright
- High brightness => triangle are too wide: a wide triangle collects a large amount of light and focuses is
  onto a signle pixel of the caustics.
- Too wide => need velocity towards the centroid of the triangle to bring its corners closer and decrease its
  surface.
- Velocity toward the centre => divergence of the velocities should be negative. That is ∇²ϕ < 0.
- The higher ∇ϕ (the actual loss given a particular target caustics), the lower ∇²ϕ (of the current estimated ϕ)
  should be.

"""
function propagate_poisson(
    ϕ::Matrix{Float64},
    ∇²ϕ::Matrix{Float64};
    clamp_correction = true,
)

    # Ensure that border conditions are as they should (= 0.0)
    fill_borders!(∇²ϕ, 0.0)
    @assert test_average_absolute(ϕ) """

        Relaxation:
            ϕ values seem too large
            max/min ∇ϕ = $(maximum(∇²ϕ)) / $(minimum(∇²ϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    @assert test_average_absolute(∇²ϕ) """

        Relaxation:
            ∇²ϕ values seem too large
            max/min ∇²ϕ = $(maximum(∇²ϕ)) / $(minimum(∇²ϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    @assert any(isnan.(ϕ)) == false "Relaxation: calling with a ϕ that contains NaN!!"
    @assert any(isnan.(∇²ϕ)) == false "Relaxation: calling with a ∇²ϕ that contains NaN!!"


    ##################################################################################################################
    # From the velocity potential field ϕ, we calculate its gradient (the velocity), then its divergence.
    # This is then compared to the intensity loss


    # Laplacian
    # Lϕ represent the approximation of the Laplacian of ϕ (second-order value)
    Lϕ = laplacian(ϕ)


    # In the case of the velocity potential, δ is the difference between the divergence of the luminosity error
    # and the current one for the current field ϕ (divergence of its gradient).

    # High luminosity error ∇²ϕ means that the triangles are too bright = too large.
    # The Laplacian Lϕ has to converge towards the ∇²ϕ. Difference to calculate speed of the descent.
    δ = Lϕ - ∇²ϕ
    @assert test_average_absolute(δ) """

        Relaxation:
            difference between luminosity error and estimated Laplacian seem too large
            max/min δ = $(maximum(δ)) / $(minimum(δ))
            max/min target Lϕ = $(maximum(Lϕ)) / $(minimum(Lϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    # Default correction ratio to adjust δ is ω / 4.0
    ω_clamped = ω / 4.0

    maximum_δ = maximum(abs.(δ))

    # We limit the changes of the potential to a maximum
    max_correction_ratio = 10.0

    if 1.94 / 4.0 * maximum_δ < max_correction_ratio
        correction_ratio = 1.94 / 4.0
    else
        correction_ratio = max_correction_ratio / maximum_δ
    end

    # if δ > 0, the intensity loss is too high, or the Laplacian is too high at that point.
    # => increase flow towards outside
    # => decrease the flow towards the inside.
    # => increase divergence of the flow
    # => increase the height of the potential there so that gradient are locally increasing
    # => correction has to be of the same sign as δ.
    @. ϕ[1:N_Pixel_Height, 1:N_Pixel_Width] +=
        correction_ratio * δ[1:N_Pixel_Height, 1:N_Pixel_Width]

    # # We limit the changes of the potential to a maximum
    # max_correction_max = 10.0
    # max_correction_μ = 20.0

    # abs_δ = abs.(δ)
    # maximum_δ = maximum(abs_δ)

    # # A correction factor is calculated to limit too wide changes on the field ϕ.
    # if clamp_correction
    #     μ_δ = average(abs_δ)

    #     if ω_clamped * maximum_δ > max_correction_max
    #         ω_clamped = min(ω_clamped, max_correction_max / maximum_δ)
    #         println(
    #             "Propagation clamped for maximum_δ.  $(maximum_δ) > $(max_correction_max)",
    #         )
    #     elseif ω_clamped *  μ_δ > max_correction_μ
    #         ω_clamped = min(ω_clamped, max_correction_μ / μ_δ)
    #         println("Propagation clamped for average δ.  $(μ_δ) > $(max_correction_μ)")
    #     end
    # end

    # In the case of the velocity potential,
    # if δ > 0, the Laplacian is too high at that point (∇²ϕ is too low).

    # To decrease the Laplacian ∇²ϕ, note that ∇² is a divergence (of a gradient) that represents net flow going outside.
    #
    # Decrease the Laplacian
    # => decrease a divergence
    # => decrease flow towards outside
    #
    # The flow is a gradient of the potential ϕ.
    # => decrease flow towards outside
    # => either decrease the values of the potential around that point, or increase the value of the potential at that point.
    # => Increase the value of ϕ at that point means a positive change when δ > 0.
    # δ > 0 => correction of the SAME sign.
    @. ϕ[1:N_Pixel_Height, 1:N_Pixel_Width] +=
        ω_clamped * δ[1:N_Pixel_Height, 1:N_Pixel_Width]
    @assert !any(isnan.(ϕ)) """

        Relaxation:
            Updated ϕ contains NaN!!

        """

    max_update = ω_clamped * maximum_δ
    return ϕ, max_update
end



"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh!(mesh::FaceMesh)

    println("""
    March mesh:
        max/min ϕ = $(maximum(abs.(mesh.corners.ϕ))) / $(minimum(abs.(mesh.corners.ϕ)))
        size ϕ = $(size(mesh.corners.ϕ))
        average ϕ = $(average_absolute(mesh.corners.ϕ))
        """)

    # Calculate the gradient of the velocity potential.
    # _POSTS_SIZED_
    ∇ϕᵤ, ∇ϕᵥ = ∇(mesh.corners.ϕ)

    # # Normalize the vectors
    # vector_length = zeros(Float64, N_Pixel_Height + 1, N_Pixel_Width + 1)
    # @. vector_length = sqrt(∇ϕᵤ^2 + ∇ϕᵥ^2)
    # ∇ϕᵤ[vector_length .!= 0.0] ./= vector_length[vector_length .!= 0.0]
    # ∇ϕᵥ[vector_length .!= 0.0] ./= vector_length[vector_length .!= 0.0]

    @assert !any(isnan.(∇ϕᵤ)) """

        march_mesh!:
            ∇ϕᵤ contains NaN!!

        """
    @assert !any(isnan.(∇ϕᵥ)) """

        march_mesh!:
            ∇ϕᵥ contains NaN!!

        """

    println("""
    march_mesh!:
        max/min ∇ϕᵤ = $(maximum(abs.(∇ϕᵤ))) / $(minimum(abs.(∇ϕᵤ)))
        size ∇ϕᵤ = $(size(∇ϕᵤ))
        average ∇ϕᵤ = $(average_absolute(∇ϕᵤ))
        max/min ∇ϕᵥ = $(maximum(abs.(∇ϕᵥ))) / $(minimum(abs.(∇ϕᵥ)))
        size ∇ϕᵥ = $(size(∇ϕᵥ))
        average ∇ϕᵥ = $(average_absolute(∇ϕᵥ))
        """)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.
    # _POSTS_SIZED_
    mesh.corners.vr .= -∇ϕᵤ
    mesh.corners.vc .= -∇ϕᵥ

    # Just in case...
    reset_border_values!(mesh.corners)

    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)

    # # Get the time, at that velocity, for the area of the triangle to become nil.
    # # _FENCES_SIZED_
    # top_triangles = [triangle3D(mesh, row, col, :top) for row ∈ 1:N_Pixel_Height, col ∈ 1:N_Pixel_Width]
    # bot_triangles = [triangle3D(mesh, row, col, :bottom) for row ∈ 1:N_Pixel_Height, col ∈ 1:N_Pixel_Width]

    # # Get the maximum moximum move along the distance of the gradient
    # # _POSTS_SIZED_
    # max_move = zeros(Float64, N_Pixel_Height + 1, N_Pixel_Width + 1)
    # max_move[1:N_Pixel_Height, 1:N_Pixel_Width] = min.(find_maximum_t.(top_triangles), find_maximum_t.(bot_triangles))
    # @. max_move[isnan(max_move)] = 0.0


    # println("""
    # march_mesh!:
    #     max_move max/min t: $(maximum(max_move)) / $(minimum(max_move))
    #     TOP max/min t: $(maximum(find_maximum_t.(top_triangles))) / $(minimum(find_maximum_t.(top_triangles)))
    #     BOT max/min t: $(maximum(find_maximum_t.(bot_triangles))) / $(minimum(find_maximum_t.(bot_triangles)))
    # """)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the gradient
    height, width = size(mesh.corners.ϕ)
    for row ∈ 1:height-1, col ∈ 1:width-1
        list_triangles =
            [triangle3D(mesh, row, col, :top), triangle3D(mesh, row, col, :bottom)]
        list_maximum_t = [
            time for time ∈ find_maximum_t.(list_triangles) if
            !isnan(time) && time > 0.0 && time < 10.0
        ]

        if !isempty(list_maximum_t)
            min_positive_t = minimum(list_maximum_t)

            # @assert all(typeof.(t1) .== Float64) && all(typeof.(t2) .== Float64) "March mesh: Maximum times are not numerical at $(row), $(col)"

            max_shift = max(abs(∇ϕᵤ[row, col]), abs(∇ϕᵥ[row, col]))

            # No shift more than 20 pixels
            max_pixel_move = 20.0
            δ = min(min_positive_t / 2.0, max_pixel_move / max_shift)
            mesh.corners.r[row, col] =
                clamp(mesh.corners.r[row, col] + δ * ∇ϕᵤ[row, col], 1, height - 1)
            mesh.corners.c[row, col] =
                clamp(mesh.corners.c[row, col] + δ * ∇ϕᵥ[row, col], 1, width - 1)
        end
    end

    # # Maximum movement along either coordinate
    # max_pixel_move = 10.0

    # # δ = zeros(Float64, N_Pixel_Height + 1, N_Pixel_Width + 1)
    # # δ = min.(max_move / 2.0, max_pixel_move)
    # δ = max_move / 2.0
    # println("Size δ = $(size(δ))")

    # @. mesh.corners.r += δ * ∇ϕᵤ
    # @. mesh.corners.c += δ * ∇ϕᵥ

    # Reset the border at the fixed values fixed coordinates.
    reset_border_values!(mesh.corners)

    println("""
    march_mesh!:
        Average mesh changes on row = $(average_absolute(mesh.corners.r - mesh_r))
        Average mesh changes on col = $(average_absolute(mesh.corners.c - mesh_c))
        Top left rows: $(mesh.corners.r[1:5, 1:5])
        Top right rows: $(mesh.corners.r[1:5, end-5:end])
    """)

    return nothing
end
