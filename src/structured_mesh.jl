# 结构化盒体网格生成器（替代 Gmsh 的 Delaunay）。
#
# 把 [0,lx]×[0,ly]×[0,lz] 切成 nx×ny×nz 个立方体，每个立方体沿体对角线 (0,0,0)-(1,1,1)
# 用 Kuhn 三角剖分拆成 6 个四面体。相比 Gmsh 的 Delaunay，厚度方向层数完全可控——
# 薄板退磁因子对厚度分辨率很敏感，而 Delaunay 对薄板会产生 sliver 单元且厚度方向不可控。

function structured_box_mesh(nx, ny, nz; lx=1.0, ly=1.0, lz=0.1)
    Nn = (nx + 1) * (ny + 1) * (nz + 1)
    nodes = zeros(Float64, Nn, 3)
    gid(i, j, k) = (k * (ny + 1) + j) * (nx + 1) + i + 1   # 1-based
    idx = 0
    for k in 0:nz, j in 0:ny, i in 0:nx
        idx += 1
        nodes[idx, :] = [lx * i / nx, ly * j / ny, lz * k / nz]
    end
    elements = zeros(Int, 6 * nx * ny * nz, 4)
    e = 0
    for k in 0:nz-1, j in 0:ny-1, i in 0:nx-1
        v0 = gid(i, j, k);         v1 = gid(i + 1, j, k)
        v2 = gid(i + 1, j + 1, k); v3 = gid(i, j + 1, k)
        v4 = gid(i, j, k + 1);     v5 = gid(i + 1, j, k + 1)
        v6 = gid(i + 1, j + 1, k + 1); v7 = gid(i, j + 1, k + 1)
        for tet in ((v0, v1, v2, v6), (v0, v2, v3, v6), (v0, v3, v7, v6),
                    (v0, v7, v4, v6), (v0, v4, v5, v6), (v0, v5, v1, v6))
            e += 1
            elements[e, :] = collect(tet)
        end
    end
    return nodes, elements
end
