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

    marginal_change = nothing
    max_update = Inf
    counter = 0
    while (abs(max_update) > 1e-6 && counter < 10_000)
        counter += 1

        print(
            """
          ================================================================================================================
          STARTING ITERATION $(counter):
              starting ϕ field = $(field_summary(mesh.corners.ϕ))

              """,
        )

        ϕ = mesh.corners.ϕ
        mesh.corners.ϕ, mesh.corners.r, mesh.corners.c, ε, max_update =
            solve_velocity_potential!(mesh, imageBW, "it$(counter)")

        marginal_change = mesh.corners.ϕ - ϕ

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
        average_absolute(marginal_change) < 1e-2 && break
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
    # Illumination only depends on the position of the corners, not their heights. ϕ is not relevant.
    lens_pixels_area = get_lens_pixels_area(mesh)

    # Positive error => the triangle needs to shrink (less light)
    ε = Float64.(lens_pixels_area - image)
    ε .-= average(ε)

    # Save the loss image as a png
    save_plot_scalar_field!(ε, "error_" * prefix, image)

    luminosity_ratio = sum(image) / sum(lens_pixels_area)
    save(
        "./examples/img_$(prefix).png",
        Gray.(clamp.(lens_pixels_area * luminosity_ratio, 0.0, 1.0)),
    )

    # Start with a clean, flat potential field.
    ϕ = zeros(Float64, height + 1, width + 1)
    println("""
            solve_velocity_potential! before loop:
                $(field_summary(lens_pixels_area, "Pixel area"))
                $(field_summary(ε, "luminosity error"))
                """)

    count = 0
    new_update = 10_000
    old_update = 2 * new_update
    while 1e-6 < new_update < old_update &&
              1e-5 < (old_update - new_update) / old_update &&
              count < 10_000

        count += 1
        old_update = new_update

        ϕ, ∇²ϕ_est, δ, new_update = CausticsEngineering.propagate_poisson(ϕ, ε)

        count % 1_000 == 0 && println("""
                                      Iteration $(count) - max_update = $(round(new_update, sigdigits=4))
                                          $(field_summary(mesh.corners.ϕ, "ϕ"))
                                          $(field_summary(∇²ϕ_est, "∇²ϕ_est"))
                                          $(field_summary(δ, "δ"))
                                          """)
    end

    # Now we need to march the mesh row,col corner locations according to this gradient.
    coord_r, coord_c = march_mesh(mesh, ϕ)

    return ϕ, coord_r, coord_c, ε, new_update
end



"""
$(SIGNATURES)

ϕ is the height map.

`march_mesh!` flexes the mesh
"""
function march_mesh(mesh::FaceMesh, ϕ::AbstractMatrix{Float64})

    # Calculate the gradient of the velocity potential.
    ∇ϕᵤ, ∇ϕᵥ = ∇(ϕ)
    # fill_borders!(∇ϕᵤ, 0.)
    # fill_borders!(∇ϕᵥ, 0.)

    # For each point in the mesh we need to figure out its velocity
    # However all the nodes located at a border will never move
    # I.e. velocity (Vx, Vy) = (0, 0) and the square of acrylate will remain of the same size.
    # reset_border_values!(mesh.corners)
    mesh_r = copy(mesh.corners.r)
    mesh_c = copy(mesh.corners.c)

    # Get the time, at that velocity, for the area of the triangle to be nil.
    # We are only interested by positive values to only move in the direction of the (opposite) gradient
    mesh.corners.vr .= -∇ϕᵤ
    mesh.corners.vc .= -∇ϕᵥ

    height, width = size(mesh)
    list_triangles = vcat(
        [
            CausticsEngineering.triangle3D(mesh, row, col, :top) for row ∈ 1:height,
            col ∈ 1:width
        ],
        [
            CausticsEngineering.triangle3D(mesh, row, col, :bottom) for row ∈ 1:height,
            col ∈ 1:width
        ],
    )

    list_maximum_t = [
        time for time ∈ CausticsEngineering.find_maximum_t.(list_triangles) if
        !isnothing(time) && time > 0.0
    ]

    δ = minimum(list_maximum_t) / 3.0
    mesh_r .-= δ * ∇ϕᵤ
    mesh_c .-= δ * ∇ϕᵥ

    println(
        """

    March mesh with correction_ratio δ = $(δ)
        $(field_summary(∇ϕᵤ, "Vu"))
        $(field_summary(∇ϕᵥ, "Vu"))
        $(field_summary(mesh_r - mesh.corners.r , "new mesh changes on row"))
        $(field_summary(mesh_c - mesh.corners.c, "new mesh changes on col"))
        $(field_summary(mesh_r - mesh.corners.rows_numbers , "total mesh changes on row"))
        $(field_summary(mesh_c - mesh.corners.cols_numbers, "total mesh changes on col"))

        """,
    )

    return mesh_r, mesh_c
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
    # Ensure that border conditions are as they should (= 0.0)
    height, width = size(ε)

    # Laplacian
    # Lϕ represent the approximation of the Laplacian of ϕ (second-order value)
    ∇²ϕ_est = laplacian(ϕ)

    # In the case of the velocity potential, δ is the difference between the divergence of the luminosity error
    # and the current one for the current field ϕ (divergence of its gradient).

    # High luminosity error ∇²ϕ means that the triangles are too bright = too large.
    # The Laplacian Lϕ has to converge towards the ∇²ϕ. Difference to calculate speed of the descent.
    δ = ∇²ϕ_est - ε

    # Correction ratio to avoid setting the triangle to nothing
    δ *= 0.20

    # In the case of the velocity potential,
    # if δ > 0, the Laplacian is too high at that point (∇²ϕ is too low).
    ϕ_res = zeros(Float64, size(ϕ))
    ϕ_res[1:end-1, 1:end-1] = ϕ[1:end-1, 1:end-1] + δ
    ϕ_res .-= average(ϕ_res)

    return ϕ_res, ∇²ϕ_est, δ, maximum(abs.(δ))
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

    ϕ = mesh.corners.ϕ

    max_update = Inf
    old_update = Inf
    for i ∈ 1:typemax(Int64)
        old_update = max_update
        ϕ, _, _, max_update = propagate_poisson(ϕ, -divergence_direction)

        i % 1_000 == 0 &&
            println("Convergence horiz. Δ  at counter = $(i) max update = $(max_update)")
        if abs(old_update - max_update) / old_update < 1e-5 || abs(max_update) <= 1e-6
            println(
                "Convergence horiz. Δ stopped with max_update of $(max_update) at counter = $(i)",
            )
            break
        end
    end

    return ϕ, max_update
end
