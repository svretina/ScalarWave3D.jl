using ScalarWave3D
using Documenter

DocMeta.setdocmeta!(ScalarWave3D, :DocTestSetup, :(using ScalarWave3D); recursive=true)

makedocs(;
    modules=[ScalarWave3D],
    authors="Stamatis Vretinaris",
    sitename="ScalarWave3D.jl",
    format=Documenter.HTML(;
        canonical="https://svretina.github.io/ScalarWave3D.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/svretina/ScalarWave3D.jl",
    devbranch="master",
)
