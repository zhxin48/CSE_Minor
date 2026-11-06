# 薄板退磁收敛测试：用结构化网格（厚度方向 nz 层）算离散退磁因子 Nx，与 Aharoni 解析对比。
#
# 只算表面项 SS + SV（均匀 M 时 VV=0），所以比全流程快很多。
# 用法：julia --project=. src/scripts/converge_thickness.jl
# 关键：随着 nz（厚度方向层数）增加，Nx_discrete 应收敛到 Nx_Aharoni ≈ 0.0975（ratio → 1）。

include(joinpath(@__DIR__, "..", "structured_mesh.jl"))
include(joinpath(@__DIR__, "..", "volume.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "surface_stiffness.jl"))
include(joinpath(@__DIR__, "..", "solve.jl"))
using LinearAlgebra
using Printf

function discrete_Nx(nodes, elements)
    ndof = 3 * size(nodes, 1)
    Nn = size(nodes, 1)
    A = zeros(ndof, ndof)
    assemble_surface_stiffness!(A, elements, nodes; rtol=1e-3)
    u = zeros(ndof); u[1:3:end] .= 1.0
    w = A * u
    vol = zeros(Nn)
    for alpha in axes(elements, 1)
        inode = elements[alpha, :]
        x = [nodes[inode[i], :] for i in 1:4]
        vol[inode] .+= tet_volume(x...) / 4
    end
    return sum(w[1:3:end]) / sum(vol)
end

function main()
    chi = 100.0
    Nx_a = demag_factors(1.0, 1.0, 0.1)[1]          # Aharoni 解析 Nx
    println("Aharoni 解析 Nx = ", Nx_a, "  对应 mean Mx = ", chi / (1 + chi * Nx_a), "（χ=", chi, ", H0=1）\n")
    println("nx ny nz   Ne     Nx_discrete   Nx_Aharoni   ratio")
    for (nx, ny, nz) in ((6, 6, 2), (6, 6, 3), (6, 6, 4), (6, 6, 5))
        nodes, elements = structured_box_mesh(nx, ny, nz)
        Nx_d = discrete_Nx(nodes, elements)
        @printf("%d %d %d  %5d   %.5f      %.5f      %.3f\n",
                nx, ny, nz, size(elements, 1), Nx_d, Nx_a, Nx_d / Nx_a)
        flush(stdout)
    end
end
main()
