module BoundaryConditions

using LinearAlgebra
using StaticArrays

## This matrix is setup to implement absorbing/radiative BCs
@inline function Radiative(n::AbstractVector{T}) where {T}
    @fastmath @inbounds begin
        mat = @SMatrix [1/2 -n[1]/2 -n[2]/2 -n[3]/2
            -n[1]/2 n[1]^2/2+n[2]^2+n[3]^2 -n[1]*n[2]/2 -n[1]*n[3]/2
            -n[2]/2 -n[1]*n[2]/2 n[1]^2+n[3]^2+n[2]^2/2 -n[2]*n[3]/2
            -n[3]/2 -n[1]*n[3]/2 -n[2]*n[3]/2 n[1]^2+n[2]^2+n[3]^2/2]
    end
    return mat
end

@inline function Reflecting(n::AbstractVector{T}) where {T}
    @fastmath @inbounds begin
        mat = @SMatrix [0.0 0.0 0.0 0.0
            -n[1] 1.0 0.0 0.0
            -n[2] 0.0 1.0 0.0
            -n[3] 0.0 0.0 1.0]
    end
    return mat
end

@inline function apply_radiative_boundary_conditions(dU, i, j, k, N)
    n = normal_vector(i, j, k, N)
    return Radiative(n) * dU
end

@inline function apply_reflecting_boundary_conditions(dU, i, j, k, N)
    n = normal_vector(i, j, k, N)
    return Reflecting(n) * dU
end

@inline function normal_vector(i, j, k, N)
    vec = @SVector [i == 1 ? -1 : i == N ? 1 : 0,
        j == 1 ? -1 : j == N ? 1 : 0,
        k == 1 ? -1 : k == N ? 1 : 0]
    return vec / norm(vec)
end

end
