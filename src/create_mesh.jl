
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

    # imageBW is normalised to the same (sort of) _energy_ as the original image.
    imageBW /= average(imageBW)

    max_update = Inf
    counter = 0
    while (abs(max_update) > 1e-5 && counter < 5_000)
        counter += 1

        println(
            "\nSTARTING VERTICAL ITERATION $(counter) --------------------------------------------",
        )

        max_update, error_luminosity =
            solve_velocity_potential!(mesh, imageBW, "it$(counter)")
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

    # Start with a new clean field
    height, width = size(mesh)

    # Get the area of each individual pixel as stretch/shrunk on the lens. Area = energy.
    # _FENCES_SIZED_
    mesh.corners.ϕ = zeros(Float64, height + 1, width + 1)
    lens_pixels_area = get_lens_pixels_area(mesh)

    # Positive error => the triangle needs to shrink (less light)
    error_luminosity = Float64.(lens_pixels_area - image)
    error_luminosity .-= average(error_luminosity)

    # Save the loss image as a png
    save_plot_scalar_field!(error_luminosity, prefix * "_loss", image)

    println("""

            solve_velocity_potential! before loop:
                $(field_summary(lens_pixels_area, "Pixel area"))
                $(field_summary(error_luminosity, "luminosity error"))
                """)


    max_update = Inf
    old_update = Inf
    for i ∈ 1:typemax(Int64)
        old_update = max_update
        max_update = CausticsEngineering.propagate_poisson!(mesh.corners.ϕ, error_luminosity)
        i % 500 == 0 && println(
            """ solve_velocity_potential! during loop at count $(i) with max_update = $(round(max_update, sigdigits=4)):
                    $(field_summary(mesh.corners.ϕ, "mesh.corners.ϕ"))
            """,
        )

        if abs(old_update - max_update) <= 1e-6 || abs(max_update) <= 1e-4
                println(""" solve_velocity_potential! at count $(i) with
                                $(field_summary(error_luminosity, "luminosity error"))
                                $(field_summary(mesh.corners.ϕ, "mesh.corners.ϕ"))
                                max_update = $(max_update)
                                """)
                break
        end
    end

    # Now we need to march the mesh row,col corner locations according to this gradient.
    march_mesh!(mesh)
    plot_as_quiver(mesh, n_steps = 50, scale = 1.0, max_length = height / 20)

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

    # Normalised divergence
    # divergence_direction is _FENCES_SIZED_
    divergence_direction = div_row + div_col
    # divergence_direction .-= average(divergence_direction)


    max_update = Inf
    old_update = Inf
    for i ∈ 1:typemax(Int64)
        old_update = max_update
        max_update = CausticsEngineering.propagate_poisson!(mesh.corners.ϕ, error_luminosity)
        i % 100 == 0 && println("Convergence horiz. Δ max update = $(max_update)")
        (abs(old_update - max_update) <= 1e-6 || abs(max_update) <= 1e-4) && break
    end

    println("Convergence horiz. Δ stopped with max_update of $(max_update) at counter = $(counter)")

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

In order of how to think about the flow of the program for the velocity potential,
    - if δ > 0, the Laplacian is too high at that point (∇²ϕ is too low).
    - To decrease the Laplacian ∇²ϕ, note that ∇² is a divergence (of a gradient) that represents net flow going outside.

    - Decrease the Laplacian
    - => decrease a divergence
    - => decrease flow towards outside

    - The flow is a gradient of the potential ϕ.
    - => decrease flow towards outside
    - => either decrease the values of the potential around that point, or increase the value of the potential at that point.
    - => Increase the value of ϕ at that point means a positive change when δ > 0.
    - δ > 0 => correction of the SAME sign.
"""
function propagate_poisson!(ϕ::Matrix{Float64}, ∇²ϕ::Matrix{Float64})
    # Ensure that border conditions are as they should (= 0.0)
    height, width = size(∇²ϕ)
    fill_borders!(∇²ϕ, 0.0)

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
    ϕ[1:end-1, 1:end-1] .+= δ
    ϕ .+= average(δ)

    # ##################################################################################################################
    # # From the velocity potential field ϕ, we calculate its gradient (the velocity), then its divergence.
    # # This is then compared to the intensity loss
    # max_δ = 0.0
    # for row ∈ 1:height, col ∈ 1:width
    #     ∇²ϕ_est =
    #         ϕ[max(1, row - 1), col] +
    #         ϕ[min(height, row + 1), col] +
    #         ϕ[row, max(1, col - 1)] +
    #         ϕ[row, min(width, col + 1)] - 4.0 * ϕ[row, col]

    #     δ = ∇²ϕ_est - ∇²ϕ[row, col]
    #     δ /= 2.0

    #     ϕ[row, col] += δ

    #     max_δ = max(max_δ, abs(δ))
    # end

    return maximum_δ
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
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain of the same size.
    # reset_border_values!(mesh.corners)
    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the gradient
    mesh.corners.vr .= -∇ϕᵤ
    mesh.corners.vc .= -∇ϕᵥ

    height, width = size(mesh)

    list_triangles = vcat(
        [
            CausticsEngineering.triangle3D(mesh, row, col, :top) for row ∈ 1:height, col ∈ 1:width
        ],
        [
            CausticsEngineering.triangle3D(mesh, row, col, :bottom) for row ∈ 1:height, col ∈ 1:width
        ],
    )

    list_maximum_t = [time for time ∈ CausticsEngineering.find_maximum_t.(list_triangles) if !isnothing(time) && time > 0.]
    minimum(list_maximum_t)
    maximum(list_maximum_t)

    δ = minimum(list_maximum_t) / 2.0

    mesh.corners.r += δ * ∇ϕᵤ
    mesh.corners.c += δ * ∇ϕᵥ

    # for row ∈ 1:height, col ∈ 1:width
    #     list_triangles = [triangle3D(mesh, row, col, :top), triangle3D(mesh, row, col, :bottom) ]

    #     if !isempty(list_triangles)
    #         list_maximum_t = [time for time ∈ find_maximum_t.(list_triangles) if !isnothing(time)]

    #         if !isempty(list_maximum_t)
    #             δ = minimum(list_maximum_t) / 8.0

    #             mesh.corners.r[row, col] -= δ * ∇ϕᵤ[row, col]
    #             mesh.corners.c[row, col]  -= δ * ∇ϕᵥ[row, col]
    #         end
    #     end
    # end

    println("""
        March mesh:
            $(field_summary(∇ϕᵤ, "Vu"))
            $(field_summary(∇ϕᵥ, "Vu"))
            $(field_summary(mesh.corners.r - mesh_r, "mesh changes on row"))
            $(field_summary(mesh.corners.c - mesh_c, "mesh changes on col"))

            """)
    return nothing
end
