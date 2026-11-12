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

# Return the triangular faces on the top surface z = maximum(z).  A valid
# tetrahedral mesh cannot have an interior face on this outer plane, so each
# matching local face is one triangle of the visible plate surface.
function _top_surface_triangles(nodes, elements)
    zmax = maximum(nodes[:, 3])
    scale = max(maximum(nodes[:, 3]) - minimum(nodes[:, 3]), 1.0)
    atol = 100 * eps(Float64) * scale
    local_faces = ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))
    triangles = Set{NTuple{3,Int}}()
    for alpha in axes(elements, 1), face in local_faces
        ids = (elements[alpha, face[1]], elements[alpha, face[2]], elements[alpha, face[3]])
        all(abs(nodes[i, 3] - zmax) <= atol for i in ids) || continue
        push!(triangles, Tuple(sort(collect(ids))))
    end
    isempty(triangles) && error("未在 z=max 平面上找到三角形表面")
    return collect(triangles)
end

# Piecewise-linear P1 interpolation on the triangulated top surface.
function _interpolate_top_surface(nodes, triangles, nodal_values, xs, ys)
    field = fill(NaN, length(ys), length(xs))
    tol = 1e-10
    for triangle in triangles
        i, j, k = triangle
        x1, y1 = nodes[i, 1], nodes[i, 2]
        x2, y2 = nodes[j, 1], nodes[j, 2]
        x3, y3 = nodes[k, 1], nodes[k, 2]
        denominator = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
        abs(denominator) > eps(Float64) || continue

        ix1 = max(1, searchsortedfirst(xs, min(x1, x2, x3)))
        ix2 = min(length(xs), searchsortedlast(xs, max(x1, x2, x3)))
        iy1 = max(1, searchsortedfirst(ys, min(y1, y2, y3)))
        iy2 = min(length(ys), searchsortedlast(ys, max(y1, y2, y3)))
        for iy in iy1:iy2, ix in ix1:ix2
            x, y = xs[ix], ys[iy]
            lambda1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / denominator
            lambda2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / denominator
            lambda3 = 1 - lambda1 - lambda2
            if lambda1 >= -tol && lambda2 >= -tol && lambda3 >= -tol
                field[iy, ix] = lambda1*nodal_values[i] +
                                lambda2*nodal_values[j] +
                                lambda3*nodal_values[k]
            end
        end
    end
    return field
end

"""
    plot_field_components(nodes, elements, u; outdir="plots", resolution=301)

Interpolate the P1 solution on the top surface of the plate and save three
continuous top-view colour maps: `field_Mx.png`, `field_My.png`, and
`field_Mz.png`.
"""
function plot_field_components(nodes, elements, u;
                               outdir=joinpath(@__DIR__, "..", "plots"),
                               resolution=301)
    u = vec(u)
    Nn = size(nodes, 1)
    length(u) == 3Nn || error("解向量长度 $(length(u)) 与 3Nn=$(3Nn) 不一致")
    triangles = _top_surface_triangles(nodes, elements)
    xs = range(extrema(nodes[:, 1])...; length=resolution)
    ys = range(extrema(nodes[:, 2])...; length=resolution)
    components = (u[1:3:end], u[2:3:end], u[3:3:end])
    labels = ("Mx", "My", "Mz")
    filenames = ("field_Mx.png", "field_My.png", "field_Mz.png")

    mkpath(outdir)
    gr()
    for index in 1:3
        values = components[index]
        field = _interpolate_top_surface(nodes, triangles, values, xs, ys)
        finite_values = filter(isfinite, vec(field))
        isempty(finite_values) && error("表面插值没有产生有限值")
        colour_limits = if index == 1
            extrema(finite_values)
        else
            limit = maximum(abs, finite_values)
            (-limit, limit)
        end
        colour_map = index == 1 ? :viridis : :RdBu
        p = heatmap(xs, ys, field;
                    color=colour_map, clims=colour_limits,
                    xlabel="x", ylabel="y", aspect_ratio=:equal,
                    title="Top view: $(labels[index])",
                    colorbar_title=labels[index],
                    size=(760, 620), dpi=180)
        path = joinpath(outdir, filenames[index])
        savefig(p, path)
        println("分量染色图已保存：", path)
    end
    return nothing
end
