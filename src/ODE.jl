module ODE

const USE_GPU = false

using StaticArrays
using Polyester
using ParallelStencil
using ImplicitGlobalGrid
import MPI
using ..BoundaryConditions

const proj_path = pkgdir(ODE)
if occursin("cn", gethostname()) || occursin("mn", gethostname())
    const output_path = "/gpfs/svretinaris/ScalarWave/runs/"
else
    const output_path = string(proj_path, "/output/")
end

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

function base_path(folder)
    return string(output_path, folder)
end

function sims_path(folder)
    return string(output_path, folder, "/sims/")
end

@inline function is_on_boundary(i, j, k, N)
    if i == N || j == N || k == N || i == 1 || j == 1 || k == 1
        return true
    else
        return false
    end
end

@parallel_indices (i, j, k) function rhs!(du::Data.Array,
                                          U::Data.Array,
                                          params,
                                          t::Data.Number)
    h = params.h
    N = params.N
    @fastmath @inbounds begin
        h2 = 2h
        N1 = N - 1
        Π = @view U[:, :, :, 2]
        Ψx = @view U[:, :, :, 3]
        Ψy = @view U[:, :, :, 4]
        Ψz = @view U[:, :, :, 5]

        dtΦ = @view du[:, :, :, 1]
        dtΠ = @view du[:, :, :, 2]
        dtΨx = @view du[:, :, :, 3]
        dtΨy = @view du[:, :, :, 4]
        dtΨz = @view du[:, :, :, 5]

        #calculate RHS everywhere
        if i == 1
            dxΨx = (Ψx[2, j, k] - Ψx[1, j, k]) / h
            dxΠ = (Π[2, j, k] - Π[1, j, k]) / h
        elseif i == N
            dxΨx = (Ψx[end, j, k] - Ψx[N1, j, k]) / h
            dxΠ = (Π[end, j, k] - Π[N1, j, k]) / h
        else
            dxΠ = (Π[i + 1, j, k] - Π[i - 1, j, k]) / h2
            dxΨx = (Ψx[i + 1, j, k] - Ψx[i - 1, j, k]) / h2
        end
        if j == 1
            dyΨy = (Ψy[i, 2, k] - Ψy[i, 1, k]) / h
            dyΠ = (Π[i, 2, k] - Π[i, 1, k]) / h
        elseif j == N
            dyΨy = (Ψy[i, end, k] - Ψy[i, N1, k]) / h
            dyΠ = (Π[i, end, k] - Π[i, N1, k]) / h
        else
            dyΠ = (Π[i, j + 1, k] - Π[i, j - 1, k]) / h2
            dyΨy = (Ψy[i, j + 1, k] - Ψy[i, j - 1, k]) / h2
        end
        if k == 1
            dzΨz = (Ψz[i, j, 2] - Ψz[i, j, 1]) / h
            dzΠ = (Π[i, j, 2] - Π[i, j, 1]) / h
        elseif k == N
            dzΨz = (Ψz[i, j, end] - Ψz[i, j, N1]) / h
            dzΠ = (Π[i, j, end] - Π[i, j, N1]) / h
        else
            dzΠ = (Π[i, j, k + 1] - Π[i, j, k - 1]) / h2
            dzΨz = (Ψz[i, j, k + 1] - Ψz[i, j, k - 1]) / h2
        end
        dtΦ[i, j, k] = Π[i, j, k]
        dtΠ[i, j, k] = dxΨx + dyΨy + dzΨz
        dtΨx[i, j, k] = dxΠ
        dtΨy[i, j, k] = dyΠ
        dtΨz[i, j, k] = dzΠ
        # Apply boundary conditions
        if is_on_boundary(i, j, k, N)
            ijk = @SVector [dtΠ[i, j, k], dtΨx[i, j, k], dtΨy[i, j, k],
                            dtΨz[i, j, k]]
            temp = BoundaryConditions.apply_radiative_boundary_conditions(ijk, i, j,
                                                                          k, N)
            dtΠ[i, j, k] = temp[1]
            dtΨx[i, j, k] = temp[2]
            dtΨy[i, j, k] = temp[3]
            dtΨz[i, j, k] = temp[4]
        end
    end
    return nothing
end

end # end of module
