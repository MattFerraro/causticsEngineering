
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

    # imageBW is normalised to have 1 unit per pixel on average.
    imageBW /= average(imageBW)

    marginal_change = nothing
    max_update = Inf
    list_max_update = []
    counter = 0
    while (abs(max_update) > 1e-6 && counter < 5)
        counter += 1

        print(
            """
            ================================================================================================================
            STARTING ITERATION $(counter):
                starting ϕ field = $(field_summary(mesh.corners.ϕ))

            """,
        )

        ε, max_update = solve_velocity_potential!(mesh, imageBW, "it$(counter)")
        marginal_change = mesh.corners.ϕ - start_ϕ

        print(
            """

            RESULT AT ITERATION $(counter):
                Luminosity error = $(field_summary(ε))
                Vertical move max update = $(max_update)
                Marginal change = $(field_summary(marginal_change))
                end ϕ field = $(field_summary(mesh.corners.ϕ))
            ================================================================================================================

            """,
        )
        plot_as_quiver(mesh, n_steps = 50, scale = height / 10, max_length = height / 20)
        average_absolute(marginal_change) < 1e-3 && break
    end

    println("\nSTARTING HORIZONTAL DISPERSION ---")
    ϕ_next, max_update = solve_horizontal_dispersion(mesh, imageBW; f = Focal_Length)

    # Move the around a nil average.
    mesh.corners.ϕ .= ϕ_next .- average(ϕ_next)

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

    # Illumination only depends on the position of the corners, not their heights. ϕ is not relevant.
    # Get the area of each individual pixel as stretch/shrunk on the lens. Area = energy.
    # _FENCES_SIZED_
    lens_pixels_area = get_lens_pixels_area(mesh)

    # Positive error => the triangle needs to shrink (less light). Enforce nil average error.
    ε = zeros(Float64, height + 1, width + 1)
    ε[1:end-1, 1:end-1] = Float64.(lens_pixels_area - image)
    ε .-= average(ε)

    # Start with a clean, flat potential field.
    ϕ_next = zeros(Float64, height + 1, width + 1)
    count = 0
    while 1e-5 < new_update && count < 3_000
        count += 1

        ϕ_next, new_update = propagate_poisson(ϕ_next, ε)
        end
    end

    # Now we need to march the mesh row,col corner locations according to this gradient.
    r_next, c_next, δ = march_mesh(mesh, ϕ_next)
    mesh.corners.r .= r_next
    mesh.corners.c .= c_next
    mesh.corners.ϕ .= ϕ_next

    return ε, new_update
end



"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh` flexes the mesh
"""
function march_mesh(mesh, mesh_ϕ)

    # Calculate the gradient of the velocity potential and reverse its direction.
    vr, vc = ∇(mesh_ϕ)
    mesh.corners.vr .= -vr
    mesh.corners.vc .= -vc

    # Clean up the mesh borders
    # reset_border_values!(mesh.corners)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain of the same size.

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the (opposite) gradient
    height, width = size(mesh)
    list_triangles = vcat(
        [triangle3D(mesh, row, col, :top) for row ∈ 1:height, col ∈ 1:width],
        [triangle3D(mesh, row, col, :bottom) for row ∈ 1:height, col ∈ 1:width],
    )

    list_maximum_t = [t for t ∈ find_maximum_t.(list_triangles) if !isnothing(t) && t > 0.0]

    δ = minimum(list_maximum_t) / 1.5
    new_r = mesh.corners.r - δ * vr
    new_c = mesh.corners.c - δ * vc

    return new_r, new_c, δ
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
function propagate_poisson(ϕ::Matrix{Float64}, ε::Matrix{Float64})

    height, width = size(ϕ)

    ∇²ϕ_est = zeros(Float64, height, width)
    δ = zeros(Float64, height, width)
    new_ϕ = ϕ

    ω = 1.99

    for row = 1:height, col = 1:width
        val_up = row == 1 ? 0.0 : new_ϕ[row-1, col]
        val_down = row == height ? 0.0 : new_ϕ[row+1, col]
        val_left = col == 1 ? 0.0 : new_ϕ[row, col-1]
        val_right = col == width ? 0.0 : new_ϕ[row, col+1]

        ∇²ϕ_est[row, col] = val_up + val_down + val_left + val_right - 4 * new_ϕ[row, col]
        δ[row, col] = ω / 4.0 * (∇²ϕ_est[row, col] - ε[row, col])
        new_ϕ[row, col] += δ[row, col]
    end

    return new_ϕ, maximum(abs.(δ))
end


"""
$(SIGNATURES)
"""
function solve_horizontal_dispersion(mesh::FaceMesh, image; f = Focal_Length)
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
    true_H = zeros(Float64, height + 1, width + 1)
    true_H = -mesh.corners.ϕ .+ H

    # Normals.
    # N_row/N_col are _POSTS_SIZED_
    N_row = zeros(Float64, size(true_H))
    N_col = zeros(Float64, size(true_H))
    N_row = tan.(atan.(d_row ./ true_H) / (n₁ - n₂))
    N_col = tan.(atan.(d_col ./ true_H) / (n₁ - n₂))

    # We need to find the divergence of the Vector field described by Nx and Ny
    # div_row/div_col are _FENCES_SIZED_
    div_row = zeros(Float64, size(true_H))
    div_col = zeros(Float64, size(true_H))
    div_row[1:end-1, 1:end-1] = N_row[2:end, 1:end-1] .- N_row[1:end-1, 1:end-1]
    div_col[1:end-1, 1:end-1] = N_col[1:end-1, 2:end] .- N_col[1:end-1, 1:end-1]

    # Normalised divergence
    # divergence_direction is _FENCES_SIZED_
    divergence_direction = div_row + div_col

    ϕ_next = mesh.corners.ϕ
    max_update = Inf
    for i ∈ 1:1e6
        ϕ_next, max_update = propagate_poisson(ϕ_next, -divergence_direction)

        i % 250 == 1 &&
            println("Convergence horiz. Δ  at counter = $(i) max update = $(max_update)")

        if abs(max_update) <= 1e-6
            println(
                "Convergence horiz. Δ stopped at counter = $(i) with max_update of $(max_update)",
            )
            break
        end
    end

    return ϕ_next, max_update
end
