"""
$(SIGNATURES)
"""
function plot_as_quiver(
    mesh::FaceMesh;
    stride = 4,
    scale = 30,
    max_length = 2,
    flipxy = false,
    reverser = false,
    reversec = false,
)

    height, width = size(mesh)
    rs = Float64[]
    cs = Float64[]
    vrs = Float64[]
    vcs = Float64[]

    ϕ = -mesh.corners.ϕ

    vrm = reverser ? mesh.corners.vr .* scale : -mesh.corners.vr .* scale
    vcm = reversec ? -mesh.corners.vc .* scale : mesh.corners.vc .* scale

    vrm = clamp.(vrm, -max_length, max_length)
    vcm = clamp.(vcm, -max_length, max_length)

    for row = 1:stride:height, col = 1:stride:width
        reverser ? push!(rs, row) : push!(rs, -row)
        reversec ? push!(cs, -col) : push!(cs, col)

        push!(vrs, vrm[row, col])
        push!(vcs, vcm[row, col])
    end

    q =
        flipxy ? quiver(cs, rs, quiver = (vcs, vrs), aspect_ratio = :equal) :
        quiver(rs, cs, quiver = (vrs, vcs), aspect_ratio = :equal)

    display(q)
end


"""
$(SIGNATURES)
"""
function plot_as_quiver(
    ϕ;
    stride = 4,
    scale = 300,
    max_length = 2,
    flipxy = false,
    reversey = false,
    reversex = false,
)

    h, w = size(ϕ)
    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]

    for x = 1:stride:w, y = 1:stride:h
        reversex ? push!(xs, x) : push!(xs, -x)
        reversey ? push!(ys, -y) : push!(ys, y)

        p1 = ϕ[y, x]
        u = (ϕ[y, x+1] - ϕ[y, x]) * scale
        v = (ϕ[y+1, x] - ϕ[y, x]) * scale

        u = -u

        reversey && (v = -v)
        reversex && (u = -u)

        # println(u, v)
        u >= 0 ? push!(us, min(u, max_length)) : push!(us, max(u, -max_length))
        v >= 0 ? push!(vs, min(v, max_length)) : push!(vs, max(v, -max_length))
    end

    q =
        flipxy ? quiver(ys, xs, quiver = (vs, us), aspect_ratio = :equal) :
        quiver(xs, ys, quiver = (us, vs), aspect_ratio = :equal)

    display(q)
end


"""
$(SIGNATURES)
"""
function plot_velocities_as_quiver(vx, vy; stride = 4, scale = 300, max_length = 2)
    h, w = size(vx)

    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]

    for x = 1:stride:w, y = 1:stride:h
        push!(xs, x)
        push!(ys, h - y)

        u = max(vx[x, y], 0.001)
        v = max(vy[x, y], 0.001)

        push!(us, u)
        push!(vs, v)
        # println(u, ": ", v)
    end

    # readline()
    q = quiver(xs, ys, quiver = (us, vs), aspect_ratio = :equal)
    display(q)
    readline()
end



"""
$(SIGNATURES)
"""
function plot_scalar_field!(field, filename, img)
    blue = zeros(size(field))
    red = zeros(size(field))
    green = zeros(size(field))

    blue[field.>0] = field[field.>0]
    red[field.<0] = -field[field.<0]

    rgbImg = RGB.(red, green, blue)'
    save("./examples/$(filename).png", map(clamp01nan, rgbImg))
end


"""
$(SIGNATURES)
"""
function plot_scalar_field(scalar_field, suffix, img)
    normalised_D_max = scalar_field ./ maximum(scalar_field)
    normalised_D_min = scalar_field ./ minimum(scalar_field)

    blue = zeros(size(scalar_field))
    red = zeros(size(scalar_field))
    green = zeros(size(scalar_field))

    blue[scalar_field.>0] = normalised_D_max[scalar_field.>0]
    red[scalar_field.<0] = -normalised_D_min[scalar_field.<0]

    rgbImg = RGB.(red, green, blue)'
    save("./examples/$(suffix)_loss.png", map(clamp01nan, rgbImg))

    println("Saving output image:")
    println(typeof(img))
    E = Gray.(scalar_field)
    println(typeof(E))
    outputImg = img - E
    save("./examples/$(suffix)_actual.png", outputImg)
end
