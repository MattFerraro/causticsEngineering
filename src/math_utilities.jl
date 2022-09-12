
"""
$(SIGNATURES)
"""
average(m::AbstractMatrix) = (sum(m) / length(m))


"""
$(SIGNATURES)
"""
average_absolute(m::AbstractMatrix) = (sum(abs.(m)) / length(m))


"""
$(SIGNATURES)
"""
test_average_absolute(m::AbstractMatrix) = (average_absolute(m)) < 1e6


"""
$(SIGNATURES)
"""
function smallest_positive(x1::Float64, x2::Float64)
    (x1 > 0.0 && x2 < 0.0) && return x1
    (x1 < 0.0 && x2 > 0.0) && return x2
    (x1 > 0.0 && x2 > 0.0) && return min(x1, x2)
    return nothing
end


"""
$(SIGNATURES)

Approximates the gradient of a scalar field.
"""
function ∇(ϕ::AbstractMatrix{Float64})
    height, width = size(ϕ)

    # x, y only represents the first and the second variables. No reference to graphical representations.
    ∇ϕx = zeros(Float64, height, width)   # divergence on the right edge will be filled with zeros
    ∇ϕy = zeros(Float64, height, width)   # divergence on bottom edge will be filled with zeros

    ∇ϕx[begin:end-1, :] = ϕ[begin+1:end, :] - ϕ[begin:end-1, :]
    ∇ϕy[:, begin:end-1] = ϕ[:, begin+1:end] - ϕ[:, begin:end-1]

    fill_borders!(∇ϕx, 0.0)
    fill_borders!(∇ϕy, 0.0)

    return ∇ϕx, ∇ϕy
end


"""
$(SIGNATURES)

Calculate the second-order Laplace operator by convolution of a kernel.
"""
function laplacian(ϕ::AbstractMatrix{Float64})
    height, width = size(ϕ)

    # Embed matrix within a larger matrix for better vectorization and avoid duplicated code
    # Padded matrix adds 1 row/col after the size of ϕ.
    # ϕ is inserted in padded matrix within  1:1+height+1 x 1:1+width+1.
    # The rest of the padded matrix (its borders) are set at 0.0.
    padded_ϕ = zeros(Float64, 1 + height + 1, 1 + width + 1)
    padded_ϕ[begin+1:end-1, begin+1:end-1] .= ϕ[1:end, 1:end]

    # Convolution = up + down + left + right
    ∇²ϕ =
        padded_ϕ[begin+1-1:end-1-1, begin+1+0:end-1+0] +
        padded_ϕ[begin+1+1:end-1+1, begin+1+0:end-1+0] +
        padded_ϕ[begin+1+0:end-1+0, begin+1-1:end-1-1] +
        padded_ϕ[begin+1+0:end-1+0, begin+1+1:end-1+1] -
        4.0 * padded_ϕ[begin+1:end-1, begin+1:end-1]

    return ∇²ϕ[1:end-1, 1:end-1]
end
