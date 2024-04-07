module Run

using ..ODE
using ..InitialData
using ..Grids
using StaticArrays

function cartesian3D(L::Real, ncells::Int, tf::Real, cfl::Real;
                     save_every::Int=1, res::String="")
    T = Float64
    dom = [-L, L]
    grid = UniformGrid(dom, ncells)
    N = grid.ncells + 1

    domt = [zero(tf), tf]
    nt = ceil(Int64, tf / (cfl * grid.spacing))
    t = Grids.UniformGrid(domt, nt)
    # Initial Data
    params_gaussian = SVector{2,Float64}(1, 0.3)
    statevector = Array{T,4}(undef, N, N, N, 5)
    # Initialize arrays
    sample_statevector!(statevector, grid, params_gaussian)

    tspan = Tuple(T.(domt.domain))
    params = (h=spacing(spacing), N=N)

    Integrator.solve(ODE.rhs!, statevector, params, t, g, save_every, folder)
    return nothing
end

function sample_statevector!(statevector::AbstractArray{T,4},
                             grid::AbstractVector{T},
                             paramsID::SVector{4,<:Integer},
                             t::T) where {T<:Real}
    f(x1, y1, z1) = InitialData.Φ(t, x1, y1, z1, paramsID)
    ppi(x1, y1, z1) = InitialData.Π(t, x1, y1, z1, paramsID)
    psix(x1, y1, z1) = InitialData.Ψx(t, x1, y1, z1, paramsID)
    psiy(x1, y1, z1) = InitialData.Ψy(t, x1, y1, z1, paramsID)
    psiz(x1, y1, z1) = InitialData.Ψz(t, x1, y1, z1, paramsID)
    N = length(grid)
    Φ = @view statevector[:, :, :, 1]
    Π = @view statevector[:, :, :, 2]
    Ψx = @view statevector[:, :, :, 3]
    Ψy = @view statevector[:, :, :, 4]
    Ψz = @view statevector[:, :, :, 5]

    @fastmath @inbounds for k in 1:N, j in 1:N, i in 1:N
        Φ[i, j, k] = f(grid[i], grid[j], grid[k])
        Π[i, j, k] = ppi(grid[i], grid[j], grid[k])
        Ψx[i, j, k] = psix(grid[i], grid[j], grid[k])
        Ψy[i, j, k] = psiy(grid[i], grid[j], grid[k])
        Ψz[i, j, k] = psiz(grid[i], grid[j], grid[k])
    end
    return nothing
end

end # end of module
