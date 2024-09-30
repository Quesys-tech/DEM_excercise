using DEM_excercise
using Documenter

DocMeta.setdocmeta!(DEM_excercise, :DocTestSetup, :(using DEM_excercise); recursive = true)

makedocs(;
    modules = [DEM_excercise],
    authors = "Queue_sys <29653819+Quesys-tech@users.noreply.github.com> and contributors",
    sitename = "DEM_excercise.jl",
    format = Documenter.HTML(;
        canonical = "https://Queue_sys.github.io/DEM_excercise.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(;
    repo = "github.com/Queue_sys/DEM_excercise.jl",
    devbranch = "main"
)
