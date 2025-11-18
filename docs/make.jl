# docs/make.jl

using Pkg

# 1. Use the docs environment
Pkg.activate(@__DIR__)

# 2. Make sure the docs env "sees" the local package in the parent folder
Pkg.develop(path = "..")

# 3. Install any missing deps declared in docs/Project.toml
Pkg.instantiate()

# 4. Now we can load Documenter and the package
using Documenter
using MarineSystemsSim

makedocs(
    sitename = "MarineSystemsSim.jl",
    modules  = [MarineSystemsSim],
    remotes  = nothing,
    format   = Documenter.HTML(
        prettyurls = false,
        repolink   = nothing,
        edit_link  = nothing,
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
