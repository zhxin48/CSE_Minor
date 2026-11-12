# Verify convergence of the non-uniform solution's mean Mx. Usage:
# julia --project=. src/scripts/verify_mean_mx.jl
#
# 验证：板的非均匀解 mean Mx 的收敛性。用法：julia --project=. src/scripts/verify_mean_mx.jl
# 方法：去掉 VV（对近均匀解贡献极小，实测差 0.2%），用结构化网格扫网格尺寸看收敛趋势。
# 结果：mean Mx 收敛到 ≈12.2（不是均匀磁化近似的 9.3）；非均匀解有效退磁因子 0.072 < 磁测 0.0975。
include(joinpath(@__DIR__, "..", "structured_mesh.jl"))
include(joinpath(@__DIR__, "..", "volume.jl"))
include(joinpath(@__DIR__, "..", "mass.jl"))
include(joinpath(@__DIR__, "..", "load.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "surface_stiffness.jl"))
include(joinpath(@__DIR__, "..", "solve.jl"))
using LinearAlgebra, Printf

# Return the mean Mx of the self-consistent solution (optionally omitting the VV term).
function mean_Mx(nodes, elements; with_vv=true)
    ndof = 3 * size(nodes, 1)
    Nn = size(nodes, 1)
    A1 = zeros(ndof, ndof); assemble_mass!(A1, elements, nodes, 100.0)
    A2 = zeros(ndof, ndof)
    if with_vv
        assemble_stiffness!(A2, elements, nodes; precision=:fast)
    end
    assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-3)
    b = zeros(ndof); assemble_load!(b, elements, nodes, 1.0)
    u = (A1 + A2) \ b
    vol = zeros(Nn)
    for al in axes(elements, 1)
        inode = elements[al, :]; x = [nodes[inode[i], :] for i in 1:4]
        vol[inode] .+= tet_volume(x...) / 4
    end
    return sum(u[1:3:end] .* vol) / sum(vol)
end

# Check convergence of mean Mx under mesh refinement.
function main()
    # 先验证 6x6x2 下 VV 是否可忽略
    nodes, elements = structured_box_mesh(6, 6, 2)
    m_full = mean_Mx(nodes, elements; with_vv=true)
    m_novv = mean_Mx(nodes, elements; with_vv=false)
    println("6x6x2: 含 VV = ", round(m_full, digits=4), "  不含 VV = ", round(m_novv, digits=4), "\n")

    println("nx ny nz   Ne     mean Mx (无 VV)")
    for (nx, ny, nz) in ((6, 6, 2), (8, 8, 2), (10, 10, 2), (12, 12, 3))
        nodes, elements = structured_box_mesh(nx, ny, nz)
        m = mean_Mx(nodes, elements; with_vv=false)
        @printf("%d %d %d  %5d   %.4f\n", nx, ny, nz, size(elements, 1), m)
        flush(stdout)
    end
end
main()
