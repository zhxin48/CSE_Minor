# Global stiffness-matrix assembly on the full mesh + structural checks.
# Reuses assemble_stiffness! (precision=:fast = self α=β → method B, else A).
#
# 完整网格（262 四面体）全局刚度矩阵组装 + 结构验证。
# 复用 assemble_stiffness!（precision=:fast = 自作用 α=β → 方法 B，其余 → 方法 A）。
include(joinpath(@__DIR__, "..", "mesh.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
using LinearAlgebra

# Assemble A2 on the standard mesh and report symmetry, row/column sums, and statistics.
function main()
    nodes, elements = load_mesh()
    Nn = size(nodes, 1)
    Ne = size(elements, 1)
    ndof = 3 * Nn
    println("Nn=$Nn  Ne=$Ne  ndof=$ndof  总单元对=$(Ne^2)")

    A2 = zeros(ndof, ndof)

    println("开始组装（precision=:fast）...")
    t0 = time()
    assemble_stiffness!(A2, elements, nodes; precision=:fast)
    t = time() - t0
    println("组装耗时 = ", round(t, digits=1), " s")

    # --- 结构验证 ---
    sym_err = maximum(abs.(A2 - A2'))
    scale = maximum(abs.(A2))
    rowsum = maximum(abs.(sum(A2, dims=2)))
    colsum = maximum(abs.(sum(A2, dims=1)))
    nz = count(!iszero, A2)
    tr = sum(diag(A2))
    dmin, dmax = extrema(diag(A2))
    total = sum(A2)

    println("\n===== 验证指标 =====")
    println("对称性 max|A2-A2'| = ", sym_err)
    println("元素尺度 max|A2|   = ", scale)
    println("相对对称偏差       = ", sym_err / scale)
    println("行和 max           = ", rowsum)
    println("列和 max           = ", colsum)
    println("总 sum(A2)         = ", total)
    println("非零元素数         = ", nz, " / ", ndof^2, "  (稠密度 ", round(nz/ndof^2, digits=4), ")")
    println("对角和 trace       = ", tr)
    println("对角 min / max     = ", dmin, " / ", dmax)

    # 保存摘要
    open(joinpath(@__DIR__, "..", "..", "data", "global_assembly_summary.txt"), "w") do io
        println(io, "Nn=", Nn, " Ne=", Ne, " ndof=", ndof, " pairs=", Ne^2)
        println(io, "precision=fast (method=hybrid, rtol=1e-3, atol=1e-6)")
        println(io, "time_s=", t)
        println(io, "sym_err=", sym_err)
        println(io, "scale=", scale)
        println(io, "rel_sym=", sym_err/scale)
        println(io, "rowsum_max=", rowsum)
        println(io, "colsum_max=", colsum)
        println(io, "sum=", total)
        println(io, "nonzero=", nz, "/", ndof^2)
        println(io, "trace=", tr)
        println(io, "diag_min=", dmin, " diag_max=", dmax)
    end
    println("\n摘要已写入 data/global_assembly_summary.txt")
end

main()
