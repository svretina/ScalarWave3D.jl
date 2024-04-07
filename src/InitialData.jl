module InitialData

# 3D wave
@inline function Φ(t::Real, x::Real, y::Real, z::Real, params)::Real
    A = params[1]
    n = params[3]
    L = params[4]
    @fastmath nl = n / L
    return @fastmath A * sinpi(nl * x) * sinpi(nl * y) * sinpi(nl * z) * cospi(√3 * nl * t)
end

@inline function Π(t::Real, x::Real, y::Real, z::Real, params)::Real
    A = params[1]
    n = params[3]
    L = params[4]
    @fastmath nl = n / L
    @fastmath ω = √3 * nl * π
    return @fastmath -A * ω * sinpi(nl * x) * sinpi(nl * y) * sinpi(nl * z) *
                     sinpi(√3 * nl * t)
end

@inline function Ψx(t::Real, x::Real, y::Real, z::Real, params)::Real
    A = params[1]
    n = params[3]
    L = params[4]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * x) * sinpi(nl * y) * sinpi(nl * z) *
                     cospi(√3 * nl * t)
end

@inline function Ψy(t::Real, x::Real, y::Real, z::Real, params)::Real
    A = params[1]
    n = params[3]
    L = params[4]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * y) * sinpi(nl * x) * sinpi(nl * z) *
                     cospi(√3 * nl * t)
end

@inline function Ψz(t::Real, x::Real, y::Real, z::Real, params)::Real
    A = params[1]
    n = params[3]
    L = params[4]
    @fastmath nl = n / L
    return @fastmath A * nl * π * cospi(nl * z) * sinpi(nl * y) * sinpi(nl * x) *
                     cospi(√3 * nl * t)
end


end
