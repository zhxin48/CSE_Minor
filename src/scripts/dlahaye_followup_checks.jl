using DelimitedFiles
using HCubature
using InteractiveUtils
using LinearAlgebra
using Printf
using Plots

include(joinpath(@__DIR__, "..", "mesh.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "surface_stiffness.jl"))

root = normpath(joinpath(@__DIR__, "..", ".."))
outdir = joinpath(root, "data", "dlahaye_followup")
plotdir = joinpath(root, "plots")
mkpath(outdir)
mkpath(plotdir)

ref_tet = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
           [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
shift_tet(t) = tuple((ref_tet[i] .+ t for i in 1:4)...)

function warntype_text(f, argtypes)
    io = IOBuffer()
    code_warntype(io, f, argtypes; debuginfo=:none, optimize=true)
    return String(take!(io))
end

function write_warntype(name, f, argtypes)
    txt = warntype_text(f, argtypes)
    write(joinpath(outdir, "warntype_$(name).txt"), txt)
    any_lines = count(line -> occursin(r"\b(?:Any|ANY)\b", line), split(txt, '\n'))
    return any_lines
end

println("[1/4] Type-stability audit")
vtype = Vector{Float64}
mtype = Matrix{Int}
tettype = NTuple{4,Vector{Float64}}
facetype = NTuple{3,Vector{Float64}}

probe_num(alpha::tettype, beta::tettype) =
    tet_tet_integral_num(alpha, beta; rtol=1e-3, atol=1e-6, maxevals=10^6)
probe_facets(alpha::tettype, beta::tettype) =
    tet_tet_integral_facets(alpha, beta; rtol=1e-3, atol=1e-6, maxevals=10^6)
probe_self_face(face::facetype) =
    self_facet_matrix(face; rtol=1e-3, atol=1e-6, maxevals=10^6)

type_rows = [
    ("basis_grads", write_warntype("basis_grads", basis_grads,
                                    Tuple{vtype, vtype, vtype, vtype})),
    ("tet_facets", write_warntype("tet_facets", tet_facets,
                                   Tuple{vtype, vtype, vtype, vtype})),
    ("boundary_faces", write_warntype("boundary_faces", boundary_faces,
                                       Tuple{mtype, Matrix{Float64}})),
    ("tet_tet_integral_num", write_warntype("tet_tet_integral_num", probe_num,
                                             Tuple{tettype, tettype})),
    ("tet_tet_integral_facets", write_warntype("tet_tet_integral_facets", probe_facets,
                                                Tuple{tettype, tettype})),
    ("self_facet_matrix", write_warntype("self_facet_matrix", probe_self_face,
                                          Tuple{facetype})),
]
open(joinpath(outdir, "warntype_summary.csv"), "w") do io
    println(io, "function,lines_containing_Any")
    for row in type_rows
        println(io, join(row, ','))
    end
end

# Count the outer evaluations, the number of inner solves (one per outer
# evaluation), and the total number of inner-integrand evaluations.
function integrate_tet_counted(f, x1, x2, x3, x4;
                               rtol=1e-3, atol=1e-6, maxevals=10^6)
    J = hcat(x2 - x1, x3 - x1, x4 - x1)
    detJ = abs(det(J))
    integrand = function (uvw)
        u, v, w = uvw
        jac = (1 - u)^2 * (1 - v)
        xi = [u, v * (1 - u), w * (1 - u) * (1 - v)]
        return f(x1 + J * xi) * detJ * jac
    end
    return hcubature_count(integrand, zeros(3), ones(3);
                            rtol=rtol, atol=atol, maxevals=maxevals)
end

function nested_tet_counted(alpha, beta;
                            outer_rtol=1e-3, inner_rtol=1e-3,
                            atol=1e-6, maxevals=10^6)
    inner_calls = Ref(0)
    inner_evals = Ref(0)
    max_inner_evals = Ref(0)
    inner = function (rd)
        inner_calls[] += 1
        val, _, n = integrate_tet_counted(rs -> 1 / norm(rd - rs), beta...;
                                           rtol=inner_rtol, atol=atol,
                                           maxevals=maxevals)
        inner_evals[] += n
        max_inner_evals[] = max(max_inner_evals[], n)
        return val
    end
    value, estimate, outer_evals = integrate_tet_counted(inner, alpha...;
                                                          rtol=outer_rtol,
                                                          atol=atol,
                                                          maxevals=maxevals)
    return (; value, estimate, outer_evals,
            inner_calls=inner_calls[], inner_evals=inner_evals[],
            max_inner_evals=max_inner_evals[])
end

println("[2/4] Nested HCubature tolerance and evaluation counts")
count_cases = [("shared_vertex", [1.0, 0.0, 0.0]),
               ("far", [5.0, 0.0, 0.0])]
tolerance_pairs = [(1e-3, 1e-3), (1e-3, 1e-4), (1e-4, 1e-3)]
open(joinpath(outdir, "nested_hcubature_counts.csv"), "w") do io
    println(io, "case,outer_rtol,inner_rtol,value,outer_error_estimate,outer_evals,inner_calls,total_inner_evals,max_inner_evals")
    for (name, shift) in count_cases
        beta = shift_tet(shift)
        for (outer_rtol, inner_rtol) in tolerance_pairs
            result = nested_tet_counted(ref_tet, beta;
                                        outer_rtol=outer_rtol,
                                        inner_rtol=inner_rtol)
            println(io, join((name, outer_rtol, inner_rtol,
                              result.value, result.estimate,
                              result.outer_evals, result.inner_calls,
                              result.inner_evals, result.max_inner_evals), ','))
            @printf("  %-13s outer=%1.0e inner=%1.0e: outer=%d, inner total=%d\n",
                    name, outer_rtol, inner_rtol,
                    result.outer_evals, result.inner_evals)
        end
    end
end

println("[3/4] Boundary extraction checks")
nodes, elements = load_mesh(filename=joinpath(root, "data", "my_box.msh"))
bfaces = boundary_faces(elements, nodes)
face_keys = [Tuple(sort([i, j, k])) for (i, j, k, _) in bfaces]
unique_keys = length(unique(face_keys))
unit_normal_error = maximum(abs(norm(n) - 1) for (_, _, _, n) in bfaces)

outward_dots = Float64[]
for (i, j, k, n) in bfaces
    owners = [a for a in axes(elements, 1)
              if i in elements[a, :] && j in elements[a, :] && k in elements[a, :]]
    @assert length(owners) == 1
    fourth = first(setdiff(elements[only(owners), :], [i, j, k]))
    push!(outward_dots, dot(n, nodes[i, :] - nodes[fourth, :]))
end
open(joinpath(outdir, "boundary_check.txt"), "w") do io
    println(io, "nodes=$(size(nodes, 1))")
    println(io, "tetrahedra=$(size(elements, 1))")
    println(io, "boundary_faces=$(length(bfaces))")
    println(io, "unique_boundary_keys=$unique_keys")
    println(io, "max_unit_normal_error=$unit_normal_error")
    println(io, "minimum_outward_dot=$(minimum(outward_dots))")
    println(io, "all_faces_have_one_owner=true")
end

println("[4/4] Separate Mx, My, Mz panels")
u = vec(readdlm(joinpath(root, "data", "u.txt"), Float64))
nn = size(nodes, 1)
@assert length(u) == 3nn
components = (u[1:3:end], u[2:3:end], u[3:3:end])
groups = Dict{Tuple{Float64,Float64},Vector{Int}}()
for i in 1:nn
    key = (round(nodes[i, 1]; digits=12), round(nodes[i, 2]; digits=12))
    push!(get!(groups, key, Int[]), i)
end
xykeys = sort(collect(keys(groups)); by=k -> (k[2], k[1]))
x = first.(xykeys)
y = last.(xykeys)
averaged = [[sum(comp[i] for i in groups[k]) / length(groups[k]) for k in xykeys]
            for comp in components]
mag = sqrt.(averaged[1].^2 .+ averaged[2].^2 .+ averaged[3].^2)

gr()
titles = ["M_x", "M_y", "M_z", "|M|"]
values = [averaged[1], averaged[2], averaged[3], mag]
panels = [scatter(x, y; marker_z=values[i], color=:viridis, markersize=5,
                  markerstrokewidth=0.4, markerstrokecolor=:black,
                  aspect_ratio=:equal, xlabel="x", ylabel="y",
                  title=titles[i], legend=false, colorbar=true)
          for i in 1:4]
fig = plot(panels...; layout=(2, 2), size=(1050, 900), dpi=180)
savefig(fig, joinpath(plotdir, "field_components.png"))

println("Outputs written to $outdir and $(joinpath(plotdir, "field_components.png"))")
