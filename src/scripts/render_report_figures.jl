using DelimitedFiles
using LinearAlgebra
using Plots

include(joinpath(@__DIR__, "..", "mesh.jl"))

project_root = normpath(joinpath(@__DIR__, "..", ".."))
report_figures = joinpath(project_root, "The_magnetisation_of_a_metallic_object", "images")
mkpath(report_figures)

# Re-render the saved self-consistent solution without the misleading scale of
# a three-dimensional quiver plot.  Colour carries magnitude; arrows carry
# direction only and have one common visual length.
nodes, elements = load_mesh(filename=joinpath(project_root, "data", "my_box.msh"))
u = vec(readdlm(joinpath(project_root, "data", "u.txt"), Float64))
nn = size(nodes, 1)
@assert length(u) == 3nn

mx_all = u[1:3:end]
my_all = u[2:3:end]
mz_all = u[3:3:end]

# Nodes above one another in the thin direction have identical projected
# coordinates.  Average those values before plotting so that the top view has
# one arrow per (x,y) position rather than overlapping arrows from each layer.
groups = Dict{Tuple{Float64,Float64},Vector{Int}}()
for i in 1:nn
    key = (round(nodes[i, 1]; digits=12), round(nodes[i, 2]; digits=12))
    push!(get!(groups, key, Int[]), i)
end
keys_xy = sort(collect(keys(groups)); by=k -> (k[2], k[1]))
x = first.(keys_xy)
y = last.(keys_xy)
mx = [sum(mx_all[i] for i in groups[k]) / length(groups[k]) for k in keys_xy]
my = [sum(my_all[i] for i in groups[k]) / length(groups[k]) for k in keys_xy]
mz = [sum(mz_all[i] for i in groups[k]) / length(groups[k]) for k in keys_xy]
mag = sqrt.(mx .^ 2 .+ my .^ 2 .+ mz .^ 2)
inplane = sqrt.(mx .^ 2 .+ my .^ 2)
scale = 0.035
qx = scale .* mx ./ max.(inplane, eps(Float64))
qy = scale .* my ./ max.(inplane, eps(Float64))

gr()
p = scatter(x, y; marker_z=mag, color=:viridis, colorbar=true,
            markersize=5, markerstrokewidth=0.5, markerstrokecolor=:black,
            xlabel="x", ylabel="y", aspect_ratio=:equal,
            xlims=(-0.05, 1.05), ylims=(-0.05, 1.05),
            title="Self-consistent magnetisation: top view",
            colorbar_title="|M|", legend=false, size=(760, 620), dpi=180)
quiver!(p, x, y; quiver=(qx, qy),
        color=:black, linewidth=0.7, alpha=0.72)
savefig(p, joinpath(report_figures, "xin_field_top_view.png"))

# Method map for the four contributions created by the volume/surface split.
labels = ["VV", "SV", "VS", "SS"]
domains = ["Ω × Ω", "∂Ω × Ω", "Ω × ∂Ω", "∂Ω × ∂Ω"]
methods = ["Euler twice\n6D → 16 × 4D",
           "Euler once\n5D → 4 × 4D",
           "transpose of SV",
           "direct if F ≠ F′\nDuffy if F = F′"]

q = plot(; xlims=(0, 12), ylims=(0, 5.3), axis=false, grid=false,
         legend=false, size=(1100, 520), dpi=180, background_color=:white)
box!(x1, x2, y1, y2, colour) = plot!(q,
    Shape([x1, x2, x2, x1], [y1, y1, y2, y2]);
    color=colour, linecolor=:black, linewidth=1.2)

for (row, y) in enumerate((4.25, 3.15, 2.05, 0.95))
    box!(0.35, 1.75, y - 0.38, y + 0.38, :lightblue)
    box!(2.55, 5.25, y - 0.38, y + 0.38, :lightgoldenrod)
    box!(6.10, 10.85, y - 0.38, y + 0.38, :honeydew)
    annotate!(q, 1.05, y, text(labels[row], 13, :black, :center))
    annotate!(q, 3.90, y, text(domains[row], 12, :black, :center))
    annotate!(q, 8.48, y, text(methods[row], 11, :black, :center))
    plot!(q, [1.78, 2.48], [y, y]; arrow=true, color=:black, linewidth=1.2)
    plot!(q, [5.28, 6.03], [y, y]; arrow=true, color=:black, linewidth=1.2)
end
annotate!(q, 1.05, 5.00, text("term", 12, :black, :center))
annotate!(q, 3.90, 5.00, text("integration domain", 12, :black, :center))
annotate!(q, 8.48, 5.00, text("numerical treatment", 12, :black, :center))
savefig(q, joinpath(report_figures, "xin_four_term_method_map.png"))

println("Report figures written to: ", report_figures)
