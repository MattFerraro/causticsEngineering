
"""
$(SIGNATURES)

"""
function engineer_caustics(source_image)
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
    total_energy_lens = N_Pixel_Height * N_Pixel_Width * 1.0     # 1 unit of energy per pixel
    total_energy_caustics = sum(imageBW)
    average_energy_per_pixel = average(imageBW)

    # imageBW is normalised to the same (sort of) _energy_ as the original image.
    imageBW = imageBW / average_energy_per_pixel

    # start_max_update has no particular meing.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update

    error_luminosity = 1e6 * N_Pixel_Height * N_Pixel_Width
    old_error_luminosity = 2 * error_luminosity
    iteration_count = 0

    # abs((error_luminosity - old_error_luminosity) / old_error_luminosity) > 0.0 &&
    # abs((max_update - old_max_update) / old_max_update) > 0.0 &&
    while (abs(max_update) > 1e-5 && error_luminosity > 1e-5 && iteration_count < 10)
        iteration_count += 1

        println(
            "\nSTARTING VERTICAL ITERATION $(iteration_count) --------------------------------------------",
        )

        old_max_update = max_update
        old_error_luminosity = error_luminosity
        max_update, error_luminosity =
            solve_velocity_potential!(mesh, imageBW, "it$(iteration_count)")
        print(
            "Vertical move max update = $(max_update)   --  Mean error luminosity = $(error_luminosity)\n",
        )
    end

    println("\nSTARTING HORIZONTAL ITERATION ---")

    mesh.corners.ϕ, max_update = move_horizontally(mesh, imageBW; f = Focal_Length)
    println("\t Horizontal max update = $(max_update)")

    # Move the around a nil average.
    mesh.corners.ϕ .-= average(mesh.corners.ϕ)

    return mesh, imageBW
end



"""
$(SIGNATURES)

ϕ is the _velocity potential_. The velocity represents the direction towards which the mesh vertices should move.
Zones of high luminosity should be covered with more triangles and should attract more vertices.
"""
function solve_velocity_potential!(mesh, image, prefix)
    # Keep that to track how much progress is made
    mesh_ϕ = copy(mesh.corners.ϕ)

    # Get the area of each individual pixel as stretch/shrunk on the lens. Area = energy.
    # _FENCES_SIZED_
    lens_pixels_area = get_lens_pixels_area(mesh)

    # Positive error => the triangle needs to shrink (less light)
    error_luminosity = Float64.(lens_pixels_area - image)
    error_luminosity = error_luminosity .- average(error_luminosity)

    # Save the loss image as a png
    println(
        """
    Luminosity error:
        Max/min luminosity error: $(maximum(error_luminosity)) / $(minimum(error_luminosity))
        Average abs loss: $(average_absolute(error_luminosity))
        Average loss: $(average(error_luminosity))
        """,
    )
    save_plot_scalar_field!(error_luminosity, prefix * "_loss", image)


    # start_max_update has no particular meeting.
    # It is basically the initial max_update. But then is used to initialise other variables so that
    # not to stop the `while` loop on the first iteration
    start_max_update = 1_000.0
    max_update = start_max_update
    old_max_update = 2 * start_max_update
    iteration_count = 0

    while (abs(max_update) > 1e-5 && iteration_count < 1)
        iteration_count += 1
        old_max_update = max_update

        iteration_count % 10 == 0 && println(
            """

        solve_velocity_potential:
            max/min mesh.corners.ϕ = $(maximum(mesh.corners.ϕ)) / $(minimum(mesh.corners.ϕ))
            iteration count = $(iteration_count)
            max_update = $(max_update)
        """,
        )

        # Re-level ϕ around its average. It doesn't matter for potential fields
        mesh.corners.ϕ .-= average(mesh.corners.ϕ)

        # Solve the Poisson equation where the divergence of the gradient of ϕ is equal to the luminosity error.
        # High error => Too much luminosity.
        mesh.corners.ϕ, max_update = propagate_poisson(mesh.corners.ϕ, error_luminosity)
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
function move_horizontally(mesh::FaceMesh, image; f = Focal_Length)

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
    true_H = -mesh.corners.ϕ .+ H

    # Normals.
    # N_row/N_col are _POSTS_SIZED_
    N_row = zeros(Float64, size(true_H))
    N_col = zeros(Float64, size(true_H))
    N_row = tan.(atan.(d_row ./ true_H) / (n₁ - n₂))
    N_col = tan.(atan.(d_col ./ true_H) / (n₁ - n₂))

    # @. N_row[1:end, 1:end] =
    #     tan(atan(d_row[1:end, 1:end] / true_H[1:end, 1:end]) / (n₁ - n₂))
    # @. N_col[1:end, 1:end] =
    #     tan(atan(d_col[1:end, 1:end] / true_H[1:end, 1:end]) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    # div_row/div_col are _FENCES_SIZED_
    div_row = zeros(Float64, height, width)
    div_col = zeros(Float64, height, width)
    div_row[1:end, 1:end] = N_row[2:end, 1:end-1] .- N_row[1:end-1, 1:end-1]
    div_col[1:end, 1:end] = N_col[1:end-1, 2:end] .- N_col[1:end-1, 1:end-1]

    # divergence_direction = zeros(Float64, size(div_row))
    # divergence_direction is _FENCES_SIZED_
    divergence_direction = div_row + div_col

    # Normalise
    divergence_direction .-= average(divergence_direction)

    max_update = Inf
    while max_update >= 1e-5
        true_H, max_update = propagate_poisson(true_H, divergence_direction)
    end
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
function propagate_poisson(ϕ::Matrix{Float64}, ∇²ϕ::Matrix{Float64})

    # Ensure that border conditions are as they should (= 0.0)
    fill_borders!(∇²ϕ, 0.0)

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

    # Default correction ratio to adjust δ is ω / 4.0
    δ *= ω / 4.0

    # We limit the changes of the potential to a maximum
    maximum_δ = maximum(abs.(δ))

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
    @. ϕ[1:N_Pixel_Height, 1:N_Pixel_Width] += δ[1:N_Pixel_Height, 1:N_Pixel_Width]

    return ϕ, maximum_δ
end



"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh!(mesh::FaceMesh)

    # Calculate the gradient of the velocity potential.
    ∇ϕᵤ, ∇ϕᵥ = ∇(mesh.corners.ϕ)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain
    # of the same size.

    # Just in case...
    reset_border_values!(mesh.corners)

    # Basically infinity time to shrink a triangle to nothing.
    min_positive_t = Inf

    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the gradient
    mesh.corners.vr = -∇ϕᵤ
    mesh.corners.vc = -∇ϕᵥ

    list_triangles = vcat(
        [
            triangle3D(mesh, row, col, :top) for row ∈ 1:N_Pixel_Height,
            col ∈ 1:N_Pixel_Width
        ],
        [
            triangle3D(mesh, row, col, :bottom) for row ∈ 1:N_Pixel_Height,
            col ∈ 1:N_Pixel_Width
        ],
    )

    list_maximum_t = [time for time ∈ find_maximum_t.(list_triangles)]
    δ = minimum(list_maximum_t) / 2.0

    mesh.corners.r += δ * ∇ϕᵤ
    mesh.corners.c += δ * ∇ϕᵥ

    # Reset the border at the fixed values fixed coordinates.
    reset_border_values!(mesh.corners)

    println("""
        March mesh:
            max/min Vu = $(maximum(abs.(∇ϕᵤ))) / $(minimum(abs.(∇ϕᵤ)))
            max/min Vv = $(maximum(abs.(∇ϕᵥ))) / $(minimum(abs.(∇ϕᵥ)))
            Average mesh changes on row = $(average_absolute(mesh.corners.r - mesh_r))
            Average mesh changes on col = $(average_absolute(mesh.corners.c - mesh_c))

            """)


    return nothing
end
