# Render top-view colour maps of Mx, My, and Mz from the saved solution.
# Run from the project root:
#   julia --project=. src/scripts/plot_components.jl

using DelimitedFiles

include(joinpath(@__DIR__, "..", "mesh.jl"))
include(joinpath(@__DIR__, "..", "visualize.jl"))

project_root = normpath(joinpath(@__DIR__, "..", ".."))
nodes, elements = load_mesh(filename=joinpath(project_root, "data", "my_box.msh"))
u = vec(readdlm(joinpath(project_root, "data", "u.txt"), Float64))

plot_field_components(nodes, elements, u;
                      outdir=joinpath(project_root, "plots"))
