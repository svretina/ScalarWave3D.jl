module InputOutput

using WriteVTK
using ReadVTK
using YAML
using ..Grids

function get_time_grid_from_yaml(metadata_file)
    file = YAML.load_file(metadata_file)
    dom = file["Time"]["domain"]
    ncells = file["Time"]["ncells"]
    every = file["Output"]["save every"]
    time = UniformGrid(dom, div(ncells, every))
    return time
end

function writevtk(statevector, simpath, xcoord, ti, i, pvd)
    vtkpath = string(simpath, "/timestep_", i)

    vtk = vtk_grid(vtkpath, xcoord, xcoord, xcoord)

    vtk["Phi"] = statevector[:, :, :, 1]
    vtk["Pi"] = statevector[:, :, :, 2]
    vtk["Psix"] = statevector[:, :, :, 3]
    vtk["Psiy"] = statevector[:, :, :, 4]
    vtk["Psiz"] = statevector[:, :, :, 5]
    vtk["time"] = ti

    collection_add_timestep(pvd, vtk, ti)
    vtk_save(vtk)
    return nothing
end

function writevtk_initialdata(statevector, simpath, xcoord)
    pvdpath = string(simpath, "/full_simulation")
    vtkpath = string(simpath, "/timestep_0")

    pvd = paraview_collection(pvdpath)
    vtk = vtk_grid(vtkpath, xcoord, xcoord, xcoord)

    vtk["Phi"] = statevector[:, :, :, 1]
    vtk["Pi"] = statevector[:, :, :, 2]
    vtk["Psix"] = statevector[:, :, :, 3]
    vtk["Psiy"] = statevector[:, :, :, 4]
    vtk["Psiz"] = statevector[:, :, :, 5]
    vtk["time"] = 0.0

    collection_add_timestep(pvd, vtk, 0.0)
    vtk_save(vtk)
    return pvd
end

function save_pvd(pvd)
    vtk_save(pvd)
    return nothing
end

# format YAML
function write_metadata(md_path, params, save_every, simdir)
    metadata_filename = string(md_path, "/metadata.yaml")
    data = "Grid:
    domain:       [$(params.grid.domain[1]), $(params.grid.domain[2])]
    ncells:       $(params.grid.ncells)
    grid:         $(params.grid.ncells + 1)×$(params.grid.ncells + 1)×$(params.grid.ncells + 1)
    total points: $((params.grid.ncells + 1)^3)
    dx=dy=dz:     $(params.h)

Time:
    CFL:          $(params.dt/params.h)
    domain:       [0.0, $(params.ti[end])]
    ncells:       $(length(params.ti)-1)
    dt:           $(params.dt)

Output:
    output path:   $simdir
    save every:    $save_every
\n
"
    println(data)
    open(metadata_filename, "w") do file
        write(file, data)
    end
end


#     job ID:        $(params.jobid)

end # end of module
