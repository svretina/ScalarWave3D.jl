module Run

const USE_GPU = false

using ..ODE
using ..InitialData
using ..Grids
using ..Integrator
using ParallelStencil
using Polyester
using StaticArrays

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

function cartesian3D(L::Real, ncells::Int, tf::Real, cfl::Real;
                     save_every::Int=1, res::String="")
    T = Float64
    dom = [-L, L]
    grid = UniformGrid(dom, ncells)
    N = grid.ncells + 1

    domt = [zero(tf), tf]
    nt = ceil(Int64, tf / (cfl * spacing(grid)))
    t = Grids.UniformGrid(domt, nt)
    # Initial Data
    params_gaussian = SVector{3,Float64}(1, 2, L)
    # params_gaussian = SVector{3,Float64}(1, 0.2, L)

    statevector = @zeros(N, N, N, 5)
    # Initialize arrays
    sample_statevector!(statevector, coords(grid), params_gaussian, 0.0)
    @show typeof(statevector)
    params = (h=spacing(grid), N=N)

    Integrator.solve(ODE.rhs!, statevector, params, t, grid, save_every, res)
    return nothing
end

function sample_statevector!(statevector::AbstractArray{T,4},
                             grid::AbstractVector{T},
                             paramsID::SVector{3,<:Real},
                             t::T) where {T<:Real}

    # # 3D wave
    # f(x1, y1, z1) = InitialData.Φ(t, x1, y1, z1, paramsID)
    # ppi(x1, y1, z1) = InitialData.Π(t, x1, y1, z1, paramsID)
    # psix(x1, y1, z1) = InitialData.Ψx(t, x1, y1, z1, paramsID)
    # psiy(x1, y1, z1) = InitialData.Ψy(t, x1, y1, z1, paramsID)
    # psiz(x1, y1, z1) = InitialData.Ψz(t, x1, y1, z1, paramsID)

    # Gaussian
    f(x1, y1, z1) = InitialData.Gaussian(t, x1, y1, z1, paramsID)
    ppi(x1, y1, z1) = InitialData.dtGaussian(t, x1, y1, z1, paramsID)
    psix(x1, y1, z1) = InitialData.dxGaussian(t, x1, y1, z1, paramsID)
    psiy(x1, y1, z1) = InitialData.dyGaussian(t, x1, y1, z1, paramsID)
    psiz(x1, y1, z1) = InitialData.dzGaussian(t, x1, y1, z1, paramsID)

    N = length(grid)
    Φ = @view statevector[:, :, :, 1]
    Π = @view statevector[:, :, :, 2]
    Ψx = @view statevector[:, :, :, 3]
    Ψy = @view statevector[:, :, :, 4]
    Ψz = @view statevector[:, :, :, 5]

    @batch for k in 1:N
        @fastmath @inbounds for j in 1:N, i in 1:N
            Φ[i, j, k] = f(grid[i], grid[j], grid[k])
            Π[i, j, k] = ppi(grid[i], grid[j], grid[k])
            Ψx[i, j, k] = psix(grid[i], grid[j], grid[k])
            Ψy[i, j, k] = psiy(grid[i], grid[j], grid[k])
            Ψz[i, j, k] = psiz(grid[i], grid[j], grid[k])
        end
    end
    return nothing
end

end # end of module
