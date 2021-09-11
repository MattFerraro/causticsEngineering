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

    # The total energy going through the lens is equal to the amount of energy on the caustics
    total_energy_lens = height * width * 1     # 1 unit of energy per pixel
    total_energy_caustics = sum(imageBW)
    correction_ratio = sum(imageBW) / (width * height)

    # imageBW is normalised to the same (sort of) _energy_ as the original image.
    imageBW = imageBW ./ correction_ratio

    # start_max_update has no particular meing.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update
    iteration_count = 0

    while (
        1e-6 < abs(max_update) < start_max_update + 1 &&
        abs((max_update - old_max_update) / old_max_update) > 0.01 &&
        iteration_count < 1_024
    )
        iteration_count += 1

        println(
            "\nSTARTING VERTICAL ITERATION $(iteration_count) --------------------------------------------",
        )

        old_max_update = max_update
        max_update = solve_velocity_potential!(mesh, imageBW, "it$(iteration_count)")
        print("Vertical move max update = $(max_update) \n")
    end

    println("\nSTARTING HORIZONTAL ITERATION ---")
    ϕ, max_update = move_horizontally(mesh, imageBW; f = 1.0, picture_width = Caustics_Side)
    mesh.corners.ϕ .= ϕ
    println(" max update = $(max_update)")

    # Move the around a nil average.
    mean_ϕ = sum(mesh.corners.ϕ) / length(mesh.corners.ϕ)
    mesh.corners.ϕ = mesh.corners.ϕ .- mean_ϕ

    return mesh, imageBW
end



"""
$(SIGNATURES)

ϕ is the _velocity potential_. The velocity represents the direction towards which the mesh vertices should move.
Zones of high luminosity should be covered with more triangles and should attract more vertices.
"""
function solve_velocity_potential!(mesh, image, suffix)
    # Remember mesh is (will be) `grid_definition x grid_definition` just like the image
    # `grid_definition x grid_definition`, so LJ is `grid_definition x grid_definition`.

    area_distorted_corners = get_area_corners(mesh)
    @assert average_absolute(area_distorted_corners) """

           solve_velocity_potential: area_distorted_corners values seem too large
               $(maximum(area_distorted_corners)) / $(minimum(area_distorted_corners))

               """

    # Positive error => the triangle needs to shrink (less light)
    error_luminosity = Float64.(area_distorted_corners - image)
    @assert average_absolute(error_luminosity) """

            solve_velocity_potential: error_luminosity values seem too large
                max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))

                """

    # Save the loss image as a png
    println("Loss:")
    println("\tMinimum loss: $(minimum(error_luminosity))")
    println("\tMaximum loss: $(maximum(error_luminosity))")
    plot_scalar_field!(error_luminosity, suffix * "_loss", image)

    # For the purpose of Poisson, we need to divergence to increase towards the inside: negative divergence for high luminosity
    # error_luminosity = -error_luminosity
    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)
    mesh_ϕ = copy(mesh.corners.ϕ)
    height, width = size(mesh.corners.ϕ)

    # start_max_update has no particular meeting.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update
    iteration_count = 0
    new_divergence = 0.0

    while (
        1e-6 < abs(max_update) < start_max_update + 1 &&
        abs((max_update - old_max_update) / old_max_update) > 0.01 &&
        new_divergence < 100.0
    )

        # while sum(abs.(mesh.corners.ϕ - mesh_ϕ)) / (height * width) < 100
        iteration_count += 1
        iteration_count % 100 == 0 && println("Converging for intensity: $(max_update)")

        mesh.corners.ϕ[:, :] .-= sum(mesh.corners.ϕ) / (height * width)

        old_max_update = max_update
        println(
            """
        solve_velocity_potential:
            max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))
            iteration count = $(iteration_count),       max_update = $(max_update)
            new maximum divergence = $(new_divergence)
        """,
        )

        # Solve the Poisson equation where the divergence of the gradient of ϕ is equal to the luminosity error.
        # High error => Too much luminosity.
        ϕ, max_update = propagate_poisson(mesh.corners.ϕ, error_luminosity)
        mesh.corners.ϕ[:, :] .= ϕ[:, :]

        @assert average_absolute(mesh.corners.ϕ) """

                solve_velocity_potential: ϕ values seem too large
                    iteration count = $(iteration_count)
                    max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))

                    max/min error_luminosity = $(maximum(error_luminosity)) / $(minimum(error_luminosity))
                    max_update = $(max_update)

                    """
        new_divergence = maximum(
            abs.(
                mesh.corners.ϕ[begin+1:end, begin:end-1] .-
                mesh_ϕ[begin:end-1, begin:end-1] .+
                mesh.corners.ϕ[begin:end-1, begin+1:end] .-
                mesh_ϕ[begin:end-1, begin:end-1],
            ),
        )
    end


    # Relevel ϕ around its average. It doesn't matter for potential fields
    mean_ϕ = sum(mesh.corners.ϕ) / length(mesh.corners.ϕ)
    mesh.corners.ϕ .-= mean_ϕ

    println(
        """
    Convergence stopped at step $(iteration_count) with
        max_update = $(max_update)
        max/min area_distorted_corners = $(maximum(area_distorted_corners)) / $(minimum(area_distorted_corners))
        max error_luminosity = $(maximum(abs.(error_luminosity)))
        Average mesh changes on ϕ = $(sum(abs.(mesh.corners.ϕ - mesh_ϕ)) / (height * width))
        max ϕ divergence = $(maximum(abs.(mesh.corners.ϕ[begin+1:end, begin:end-1] .- mesh_ϕ[begin:end-1, begin:end-1] .+
        mesh.corners.ϕ[begin:end-1, begin+1:end] .- mesh_ϕ[begin:end-1, begin:end-1])))

        """,
    )

    # Now we need to march the x,y locations in our mesh according to this gradient!
    # Tasks are a control flow feature that allows computations to be suspended and resumed in a flexible manner.
    # We mention them here only for completeness; for a full discussion see Asynchronous Programming.
    march_mesh!(mesh)

    save_obj!(
        matrix_to_mesh(mesh.corners.ϕ * 0.01),
        "./examples/phi_$(suffix).obj",
        reverse = false,
        flipxy = false,
    )

    # plot_as_quiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reverser=false, reversec=false)
    plot_as_quiver(
        mesh,
        stride = 20,
        scale = 1.0,
        max_length = 200,
        flipxy = true,
        reverser = false,
        reversec = false,
    )

    # plot_as_quiver(ϕ * -1.0, stride=30, scale=1.0, max_length=200, flipxy=true, reversex=false, reversey=false)
    # saveObj(matrix_to_mesh(D * 10), "D_$(suffix).obj")
    save_obj!(mesh, "./examples/mesh_$(suffix).obj", flipxy = true)

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

    # The focal length is in meters. H is in corners.
    H = f / Meters_Per_Pixel

    # true_H is the real distance beteen the surface of the lens and the projection
    # H is the initial distance from the surface of the carved face to the projection screen.
    # H is constant and does not account for the change in heights due to carving.
    # Note: Higher location means closer to the caustics => negative sign
    true_H = zeros(Float64, size(mesh.corners.ϕ))
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

    # divergence_direction = zeros(Float64, size(div_row))
    divergence_direction = div_row + div_col
    @assert average_absolute(divergence_direction) """

            ALERT: horiz. Δ divergence_direction seem too large
            max/min = $(maximum(divergence_direction)) / $(minimum(divergence_direction))
            """

    true_H, max_update = propagate_poisson(true_H, -divergence_direction)

    @assert average_absolute(true_H) """

        ALERT: horiz. Δ corner_heights is too large
        max/min corner_heights = $(maximum(true_H)) / $(minimum(true_H))
        """

    println("Convergence horiz. Δ stopped with max_update of $(max_update)")

    # saveObj(matrix_to_mesh(h / 10), "./examples/heightmap.obj")
    return true_H, max_update
end


"""
$(SIGNATURES)

This function implements successive over relaxation for a matrix and its associated error matrix
There is a hardcoded assumption of Neumann boundary conditions -- that the derivative across the
boundary must be zero in all cases. See:
https://math.stackexchange.com/questions/3790299/how-to-iteratively-solve-poissons-equation-with-no-boundary-conditions

- ϕ is the potential to be solved subject to the Poisson equation and is the same size as the corners (posts)
- target is ∇ϕ and is of the size of the pixels (fences)

In order of how to think about the flow of the program:

- High ∇ϕ (intensity loss function for the velocity potential) means that the triangles are too bright
- High brightness => triangle are too wide: a wide triangle collects a large amount of light and focuses is
  onto a signle pixel of the caustics.
- Too wide => need velocity towards the centroid of the triangle to bring its corners closer and decrease its
  surface.
- Velocity toward the centre => divergence of the velocities should be negative. That is ∇²ϕ < 0.
- The higher ∇ϕ (the actual loss given a particular target caustics), the lower ∇²ϕ (of the current estimated ϕ)
  should be.

In order to reduce the error, the error informs how to change the estimated ϕ.

-
δ > 0 => divergence of ϕ needs to increase


"""
function propagate_poisson(ϕ::Matrix{Float64}, ∇²ϕ::Matrix{Float64})

    # Ensure that border conditions are as they should
    fill_borders!(∇²ϕ, 0.0)
    @assert average_absolute(ϕ) """

        Relaxation: ϕ values seem too large
            max/min ∇ϕ = $(maximum(∇²ϕ)) / $(minimum(∇²ϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    @assert average_absolute(∇²ϕ) """

        Relaxation: ∇²ϕ values seem too large
            max/min ∇²ϕ = $(maximum(∇²ϕ)) / $(minimum(∇²ϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    @assert any(isnan.(ϕ)) == false "Relaxation: calling with a ϕ that contains NaN!!"
    @assert any(isnan.(∇²ϕ)) == false "Relaxation: calling with a ∇²ϕ that contains NaN!!"

    # ϕ is of the same size as corners = number of corners + 1
    height_corners, width_corners = size(ϕ)
    height_pixels, width_pixels = size(∇²ϕ)

    # From the velocity potential field ϕ, we calculate its gradient (the velocity), then its divergence.
    # This is then compared to the intensity loss


    # Laplacian
    # Those values are all of ϕ's size and represent the  second order approximation of the
    # Laplacian of ϕ.
    Lϕ = laplacian(ϕ)
    @assert average_absolute(Lϕ) """

        Relaxation: Laplacian of ϕ a seem too large
            max/min target Lϕ = $(maximum(Lϕ)) / $(minimum(Lϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

    # δ is the difference between the divergence of the intensity loss (in the case of the velocity potential)
    # and  the current one for the current field ϕ (divergence of its gradient).
    #
    # High intensity loss ∇²ϕ means that the triangles are too bright, too large.
    # This has to converge towards the ∇²ϕ. Difference to calculate speed of the descent.
    δ = Lϕ - ∇²ϕ
    @assert average_absolute(δ) """

        Relaxation: difference between intensity loss and estimated Laplacian seem too large
            max/min δ = $(maximum(δ)) / $(minimum(δ))
            max/min target Lϕ = $(maximum(Lϕ)) / $(minimum(Lϕ))
            max/min ϕ = $(maximum(ϕ)) / $(minimum(ϕ))

            """

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
    @. ϕ[1:height_pixels, 1:width_pixels] +=
        correction_ratio * δ[1:height_pixels, 1:width_pixels]
    @assert !any(isnan.(ϕ)) "Relaxation: Updated ϕ contains NaN!!"

    max_update = correction_ratio * maximum(abs.(δ))

    return ϕ, max_update
end



"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh!(mesh::FaceMesh)

    height, width = size(mesh.corners.ϕ)

    ∇ϕᵤ, ∇ϕᵥ = ∇(mesh.corners.ϕ)
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

    # Basically infinity time to shrink a triangle to nothing.
    min_positive_t = Inf
    list_min_pos_t = zeros(Float64, (height - 2) * (width - 2))

    t1 = zeros(Float64, size(mesh))
    t2 = zeros(Float64, size(mesh))

    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)
    mesh_ϕ = copy(mesh.corners.ϕ)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the gradient
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


    # Reset the border at the fixed values fixed coordinates.
    reset_border_values!(mesh.corners)

    println(
        """
    March mesh:
        max/min Vu = $(maximum(abs.(∇ϕᵤ))) / $(minimum(abs.(∇ϕᵤ)))
        max/min Vv = $(maximum(abs.(∇ϕᵥ))) / $(minimum(abs.(∇ϕᵥ)))
        Average mesh changes on row = $(sum(abs.(mesh.corners.r - mesh_r)) / (height * width))
        Average mesh changes on col = $(sum(abs.(mesh.corners.c - mesh_c)) / (height * width))

        """,
    )

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
