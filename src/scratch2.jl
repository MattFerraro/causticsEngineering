struct P2D
    x::AbstractMatrix{Real}
    y::AbstractMatrix{Real}

    function P2D(dim1, dim2)
        x = Matrix{Float64}(undef, dim1, dim2)
        fill!(x, 0.)
        y = Matrix{Float64}(undef, dim1, dim2)
        fill!(y, 0.)

        return new(x, y)
    end
end



m = P2D(2, 3)

m.x
m.y

m.x .+= 2.0

function f1(v::P2D, a)
    v.x .+= a
end

f1(m, 3.)

