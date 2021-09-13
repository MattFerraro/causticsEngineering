
"""
$(SIGNATURES)
"""
function plot_as_quiver(
    mesh::FaceMesh;
    n_steps = 20,
    scale = 30,
    max_length = N_Pixel_Width / n_steps,
    fliprc = false,
    reverser = false,
    reversec = false,
)

    height, width = size(mesh)
    stride = Int64(ceil(max(height, width) / n_steps))

    row_s = Float64[]
    col_s = Float64[]
    v_rows = Float64[]
    v_cols = Float64[]

    # The gradient vectors reflect change of luminosity:outward flow = bigger cells = more luminosity
    # Plots are cosmetically look better when showing the light focusing to a specific location.
    # Those gradients are the opposite direction.
    ϕ = -mesh.corners.ϕ

    mat_vr = (reverser ? -1 : 1) .* mesh.corners.vr .* scale
    mat_vc = (reversec ? -1 : 1) .* mesh.corners.vc .* scale

    mat_vr = clamp.(mat_vr, -max_length, max_length)
    mat_vc = clamp.(mat_vc, -max_length, max_length)

    for row = 1:stride:height, col = 1:stride:width
        reverser ? push!(row_s, row) : push!(row_s, -row)
        reversec ? push!(col_s, -col) : push!(col_s, col)

        push!(v_rows, mat_vr[row, col])
        push!(v_cols, mat_vc[row, col])
    end

    display(
        fliprc ? quiver(row_s, col_s, quiver = (v_rows, v_cols)) :
        quiver(col_s, row_s, quiver = (v_cols, v_rows)),
    )
    return nothing
end


"""
$(SIGNATURES)
"""
function plot_as_quiver(
    ϕ;
    n_steps = 20,
    scale = 300,
    max_length = N_Pixel_Width / n_steps,
    fliprc = false,
    reverser = false,
    reversec = false,
)

    height, width = size(ϕ)
    stride = Int64(ceil(max(height, width) / n_steps))

    row_s = Float64[]
    col_s = Float64[]
    Δrow = Float64[]
    Δcol = Float64[]

    for row = 1:stride:height, col = 1:stride:width
        dr = (-(ϕ[row, col+1] - ϕ[row, col]) * scale)
        dc = (ϕ[row+1, col] - ϕ[row, col]) * scale

        reverser ? push!(row_s, col) : push!(row_s, width + 1 - col)
        reversec ? push!(col_s, height + 1 - row) : push!(col_s, row)
        reverser && (dr = -dr)
        reversec && (dc = -dc)

        push!(Δrow, clamp(dr, -maxlength, max_length))
        push!(Δcol, clamp(dc, -maxlength, max_length))
    end

    # By default, we flip to match the original image.
    display(
        fliprc ? quiver(row_s, col_s, quiver = (Δrow, Δcol)) :
        quiver(col_s, row_s, quiver = (Δcol, Δrow)),
    )

    return nothing
end


"""
$(SIGNATURES)
"""
function plot_velocities_as_quiver(
    vr,
    vc;
    n_steps = 20,
    scale = 300,
    max_length = N_Pixel_Width / n_steps,
)

    height, width = size(vr)
    stride = Int64(ceil(max(height, width) / n_steps))

    row_s = Float64[]
    col_s = Float64[]
    Δrow = Float64[]
    Δcol = Float64[]

    dr = max.(vr, 0.001)
    dc = max.(vc, 0.001)

    for row = 1:stride:height, col = 1:stride:width
        push!(row_s, height - row)
        push!(col_s, col)

        push!(Δrow, dr[row, col])
        push!(Δcol, dc[row, col])
    end

    # readline()
    display(quiver(row_s, col_s, quiver = (Δrow, Δcol), aspect_ratio = :equal))
end



"""
$(SIGNATURES)
"""
function save_plot_scalar_field!(calculated_scalar_field, filename, img)
    scalar_field = calculated_scalar_field'

    blue = zeros(size(scalar_field))
    red = zeros(size(scalar_field))
    green = zeros(size(scalar_field))

    blue[scalar_field.>0] = scalar_field[scalar_field.>0]
    red[scalar_field.<0] = -scalar_field[scalar_field.<0]

    rgbImg = RGB.(red, green, blue)'
    save("./examples/$(filename).png", map(clamp01nan, rgbImg))
end


"""
$(SIGNATURES)
"""
function save_plot_scalar_field(calculated_scalar_field, prefix, img)
    scalar_field = calculated_scalar_field'
    normalised_D_max = scalar_field ./ maximum(scalar_field)
    normalised_D_min = scalar_field ./ minimum(scalar_field)

    blue = zeros(size(scalar_field))
    red = zeros(size(scalar_field))
    green = zeros(size(scalar_field))

    blue[scalar_field.>0] = normalised_D_max[scalar_field.>0]
    red[scalar_field.<0] = -normalised_D_min[scalar_field.<0]

    rgbImg = RGB.(red, green, blue)'
    save("./examples/$(prefix)_loss.png", map(clamp01nan, rgbImg))

    println("Saving output image:")
    println(typeof(img))
    E = Gray.(scalar_field)
    println(typeof(E))
    outputImg = img - E
    save("./examples/$(prefix)_actual.png", outputImg)
end
