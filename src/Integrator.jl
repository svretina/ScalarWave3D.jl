module Integrator

const USE_GPU = false

using ..InputOutput
using ..Grids
using ..ODE
using LoopVectorization
using ParallelStencil

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end

const proj_path = pkgdir(Integrator)
if occursin("cn", gethostname()) || occursin("mn", gethostname())
    const output_path = "/gpfs/svretinaris/ScalarWave/runs/"
else
    const output_path = string(proj_path, "/output/")
end

function base_path(folder)
    return string(output_path, folder)
end

function sims_path(folder)
    return string(output_path, folder, "/sims/")
end

@inline function RK4(rhs!::F, dt::Real, reg1, reg2, reg3, reg4, params, t::Real,
                     indices::CartesianIndices, N) where {F}
    @parallel (1:N, 1:N, 1:N) rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + reg3[i] / 2
        reg2[i] = reg3[i]
    end
    @parallel (1:N, 1:N, 1:N) rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + (reg3[i] - reg2[i]) / 2
    end
    @parallel (1:N, 1:N, 1:N) rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] - reg3[i] / 2
        reg1[i] = reg1[i] + reg3[i]
        reg2[i] = reg2[i] / 6 - reg3[i]
    end
    @parallel (1:N, 1:N, 1:N) rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] + reg3[i] + reg3[i]
        reg1[i] = reg1[i] + reg2[i] + reg3[i] / 6
    end
    return nothing
end

@inline function RK4b(rhs!::F, dt::Real, reg1, reg2, reg3, reg4, params, t::Real,
                     indices::CartesianIndices, N) where {F}
    rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + reg3[i] / 2
        reg2[i] = reg3[i]
    end
    rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + (reg3[i] - reg2[i]) / 2
    end
    rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] - reg3[i] / 2
        reg1[i] = reg1[i] + reg3[i]
        reg2[i] = reg2[i] / 6 - reg3[i]
    end
    rhs!(reg4, reg1, params.h, N, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] + reg3[i] + reg3[i]
        reg1[i] = reg1[i] + reg2[i] + reg3[i] / 6
    end
    return nothing
end

function solve(rhs, statevector, params, t, grid, save_every, folder)
    # istr = get_iter_str(0, t.ncells + 1)
    base_dir = base_path(folder)
    sims_dir = sims_path(folder)
    base_str = string(sims_dir, "output_L=", Int64(grid.domain[2]), "_nc=",
                      grid.ncells)
    # dataset = string(base_str, istr, ".h5")

    if !isdir(base_dir)
        mkdir(base_dir)
    end
    if !isdir(sims_dir)
        mkdir(sims_dir)
    end

    # if isfile(dataset)
    #     rm(dataset)
    # end

    println("Output data at directory: ", sims_dir)

    # write initial data
    xcoord = coords(grid)
    h = spacing(grid)
    N = grid.ncells+1
    # pvd = InputOutput.writevtk_initialdata(statevector, sims_dir, xcoord)

    println("=========Starting time integration=========")
    print("Allocating registers for time integrator...")
    reg2 = similar(statevector)
    reg3 = similar(statevector)
    reg4 = similar(statevector)
    println("✅")

    @parallel (1:N, 1:N, 1:N) rhs(reg4, statevector, h, N, 1.2)
    @time @parallel (1:N, 1:N, 1:N) rhs(reg2, statevector, h, N, 1.2)

    ODE.rhs_batch!(reg4, statevector, h, N, 1.2)
    @time ODE.rhs_batch!(reg4, statevector, h, N, 1.2)

    println("saving every=", save_every)

    nt = t.ncells + 1
    dt = spacing(t)
    indices = CartesianIndices(statevector)
    params = (grid=grid, dt=dt, ti=coords(t), params...)
    # InputOutput.write_metadata(base_dir, params, save_every, sims_dir)
    for (i, ti) in enumerate(coords(t)[2:end])
        # println("Iteration = ", i, "/", t.ncells)
        if (i==1)  global wtime0 = Base.time()  end
        RK4(rhs, dt, statevector, reg2, reg3, reg4, params, ti, indices, N)
        # RK4(ODE.rhs_batch!, dt, statevector, reg2, reg3, reg4, params, ti, indices, N)
        # if i % save_every == 0
        #     InputOutput.writevtk(statevector, sims_dir, xcoord, ti, i, pvd)
        # end
    end
    wtime    = Base.time()-wtime0
    A_eff    = (3*2)/1e9*N*N*sizeof(Data.Number)
    wtime_it = wtime/(nt)                        # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
    println("Total steps=$nt, time=$wtime sec (@ T_eff = $(round(T_eff, sigdigits=2)) GB/s)")
    # InputOutput.save_pvd(pvd)
    return nothing
end

end # end of module
