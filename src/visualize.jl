# 可视化：mesh 线框 + 磁场分布。用 Plots（gr 后端），输出静态 PNG。
using Plots

# 提取四面体网格的所有唯一边（每个四面体 6 条边）。
function mesh_edges(elements)
    edges = Set{Tuple{Int,Int}}()
    for alpha in axes(elements, 1)
        v = elements[alpha, :]
        for (a, b) in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
            push!(edges, minmax(v[a], v[b]))
        end
    end
    return collect(edges)
end

# 把唯一边展开成带 NaN 分隔的路径坐标（Plots :path 用 NaN 断段）。
function _edge_segments(nodes, edges)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for (a, b) in edges
        append!(xs, [nodes[a, 1], nodes[b, 1], NaN])
        append!(ys, [nodes[a, 2], nodes[b, 2], NaN])
        append!(zs, [nodes[a, 3], nodes[b, 3], NaN])
    end
    return xs, ys, zs
end

"""
    plot_mesh(nodes, elements; out="plots/mesh.png")

画四面体网格 3D 线框，存 PNG。
"""
function plot_mesh(nodes, elements; out=joinpath(@__DIR__, "..", "plots", "mesh.png"))
    edges = mesh_edges(elements)
    xs, ys, zs = _edge_segments(nodes, edges)
    d = dirname(out); isempty(d) || mkpath(d)
    gr()
    plot(xs, ys, zs; seriestype=:path, legend=false, linewidth=0.6, color=:black,
         xlabel="x", ylabel="y", zlabel="z",
         title="Mesh (Nn=$(size(nodes, 1)), Ne=$(size(elements, 1)))")
    savefig(out)
    println("Mesh 图已保存：", out)
end

"""
    plot_field(nodes, elements, u; out="plots/field.png")

画磁场分布：左 3D quiver 箭头（M 矢量），右 xy 俯视散点按 |M| 着色（薄板，边缘退磁效应）。
"""
function plot_field(nodes, elements, u; out=joinpath(@__DIR__, "..", "plots", "field.png"))
    u = vec(u)
    Nn = size(nodes, 1)
    mx = [u[3i - 2] for i in 1:Nn]
    my = [u[3i - 1] for i in 1:Nn]
    mz = [u[3i]     for i in 1:Nn]
    mag = [norm([mx[i], my[i], mz[i]]) for i in 1:Nn]
    d = dirname(out); isempty(d) || mkpath(d)
    gr()
    p1 = quiver(nodes[:, 1], nodes[:, 2], nodes[:, 3],
                quiver=(mx, my, mz), legend=false,
                xlabel="x", ylabel="y", zlabel="z", title="M (quiver)")
    p2 = scatter(nodes[:, 1], nodes[:, 2], marker_z=mag, legend=false,
                 color=:viridis, xlabel="x", ylabel="y",
                 title="|M| (top view)", colorbar_title="|M|")
    p = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    savefig(p, out)
    println("磁场图已保存：", out)
end
