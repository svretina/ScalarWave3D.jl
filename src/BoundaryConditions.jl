module BoundaryConditions

using LinearAlgebra
using StaticArrays

## This matrix is setup to implement absorbing/radiative BCs
@inline function RR(n::AbstractVector)
    @fastmath @inbounds begin
        mat = @SMatrix [1/2 -n[1]/2 -n[2]/2 -n[3]/2
                        -n[1]/2 n[1]^2/2+n[2]^2+n[3]^2 -n[1] * n[2]/2 -n[1] * n[3]/2
                        -n[2]/2 -n[1] * n[2]/2 n[1]^2+n[3]^2+n[2]^2/2 -n[2] * n[3]/2
                        -n[3]/2 -n[1] * n[3]/2 -n[2] * n[3]/2 n[1]^2+n[2]^2+n[3]^2/2]
    end
    return mat
end

@inline function apply_boundary_conditions(dU, i, j, k, N)
    n = normal_vector(i, j, k, N)
    return RR(n) * dU
end

function normal_vector(i, j, k, N)
    vec = @SVector [i == 1 ? -1 : i == N ? 1 : 0,
                    j == 1 ? -1 : j == N ? 1 : 0,
                    k == 1 ? -1 : k == N ? 1 : 0]
    return vec / norm(vec)
end

end
