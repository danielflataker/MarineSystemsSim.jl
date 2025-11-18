# docs/make.jl

using Documenter
using MarineSystemsSim

DocMeta.setdocmeta!(MarineSystemsSim, :DocTestSetup, :(using MarineSystemsSim); recursive = true)

makedocs(
    sitename = "MarineSystemsSim.jl",
    modules  = [MarineSystemsSim],
    remotes  = nothing,
    format   = Documenter.HTML(
        prettyurls = false,
        repolink   = "https://github.com/danielflataker/MarineSystemsSim.jl",
        edit_link  = "main",
    ),
    pages    = [
        "Home"            => "index.md",
        "Getting started" => "getting_started.md",
        "3-DOF model"     => "model_3dof.md",
        "API reference"   => "api.md",
    ],
    checkdocs = :exports,  # only require docs for exported names
)

deploydocs(
    repo      = "github.com/danielflataker/MarineSystemsSim.jl.git",
    devbranch = "main",
    versions  = [
        "dev" => "main",   # for now, only dev docs until we start tagging
    ],
)
