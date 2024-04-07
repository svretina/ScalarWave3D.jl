module Integrator

using ..InputOutput
using LoopVectorization

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
                     indices::CartesianIndices) where {F}
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + reg3[i] / 2
        reg2[i] = reg3[i]
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i]
        reg1[i] = reg1[i] + (reg3[i] - reg2[i]) / 2
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] - reg3[i] / 2
        reg1[i] = reg1[i] + reg3[i]
        reg2[i] = reg2[i] / 6 - reg3[i]
    end
    rhs!(reg4, reg1, params, t)
    @inbounds @fastmath @avx for i in indices
        reg3[i] = dt * reg4[i] + reg3[i] + reg3[i]
        reg1[i] = reg1[i] + reg2[i] + reg3[i] / 6
    end
    return nothing
end

function solve(rhs, statevector, params, t, grid, save_every, folder)
    istr = get_iter_str(0, t.ncells + 1)
    base_dir = base_path(folder)
    sims_dir = sims_path(folder)
    base_str = string(sims_dir, "output_L=", Int64(grid.domain[2]), "_nc=",
                      grid.ncells, "_v=", v, "_")
    dataset = string(base_str, istr, ".h5")

    if !isdir(base_dir)
        mkdir(base_dir)
    end
    if !isdir(sims_dir)
        mkdir(sims_dir)
    end

    if isfile(dataset)
        rm(dataset)
    end

    println("Output data at directory: ", sims_dir)

    # write initial data
    xcoord = params.grid_coords
    pvd = InputOutput.writevtk_initialdata(statevector, sims_dir, xcoord)


    println("=========Starting time integration=========")
    print("Allocating registers for time integrator...")
    reg2 = similar(statevector)
    reg3 = similar(statevector)
    reg4 = similar(statevector)
    println("✅")

    println("saving every=", save_every)

    nt = t.ncells + 1
    dt = params.dt
    indices = CartesianIndices(statevector)

    InputOutput.write_metadata(base_dir, params, save_every, sims_dir)
    for (i, ti) in enumerate(params.ti)
        i -= 1
        if i == 0
            continue
        end
        println("Iteration = ", i, "/", t.ncells)
        @time RK4(rhs, dt, statevector, reg2, reg3, reg4, params, ti, indices)
        if i % save_every == 0
            InputOutput.writevtk(statevector, sims_dir, xcoord, ti, i, pvd)
        end
    end

    InputOutput.save_pvd(pvd)
    return nothing
end


end # end of module
