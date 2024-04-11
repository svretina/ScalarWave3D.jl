module InitialData

using ForwardDiff

# 3D wave
@inline function Φ(t::Real, x::Real, y::Real, z::Real, params)
    @inbounds A = params[1]
    @inbounds n = params[2]
    @inbounds L = params[3]
    @fastmath nl = n / L
    return @fastmath A * sinpi(nl * x) * sinpi(nl * y) * sinpi(nl * z) * cospi(√3 * nl * t)
end

@inline function Π(t::Real, x::Real, y::Real, z::Real, params)
    @inbounds A = params[1]
    @inbounds n = params[2]
    @inbounds L = params[3]
    @fastmath nl = n / L
    @fastmath ω = √3 * nl * π
    return @fastmath -A * ω * sinpi(nl * x) * sinpi(nl * y) * sinpi(nl * z) *
                     sinpi(√3 * nl * t)
end

@inline function Ψx(t::Real, x::Real, y::Real, z::Real, params)
    @inbounds A = params[1]
    @inbounds n = params[2]
    @inbounds L = params[3]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * x) * sinpi(nl * y) * sinpi(nl * z) *
                     cospi(√3 * nl * t)
end

@inline function Ψy(t::Real, x::Real, y::Real, z::Real, params)
    @inbounds A = params[1]
    @inbounds n = params[2]
    @inbounds L = params[3]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * y) * sinpi(nl * x) * sinpi(nl * z) *
                     cospi(√3 * nl * t)
end

@inline function Ψz(t::Real, x::Real, y::Real, z::Real, params)
    @inbounds A = params[1]
    @inbounds n = params[2]
    @inbounds L = params[3]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * z) * sinpi(nl * y) * sinpi(nl * x) *
                     cospi(√3 * nl * t)
end


@inline function Gaussian(t::Real, x::Real, y::Real, z::Real, params)
    A = params[1]
    σr = params[2]
    r = sqrt(x^2 + y^2 + z^2)
    if t == 0.0 && x == 0.0 && y == 0.0 && z == 0.0
        return 0.0
    else
        return A * (r - t)^3 * exp(-((r - t)^2)) / r
    end
end

@inline function dtGaussian(t::Real, x::Real, y::Real, z::Real, params)
    return ForwardDiff.derivative(t1 -> Gaussian(t1, x, y, z, params), t)
end

@inline function dxGaussian(t::Real, x::Real, y::Real, z::Real, params)
    return ForwardDiff.derivative(x1 -> Gaussian(t, x1, y, z, params), x)
end

@inline function dyGaussian(t::Real, x::Real, y::Real, z::Real, params)
    return ForwardDiff.derivative(y1 -> Gaussian(t, x, y1, z, params), y)
end

@inline function dzGaussian(t::Real, x::Real, y::Real, z::Real, params)
    return ForwardDiff.derivative(z1 -> Gaussian(t, x, y, z1, params), z)
end

end
