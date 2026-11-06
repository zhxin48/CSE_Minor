# 组装并保存三个矩阵：
#   A1 —— 质量矩阵 A^(1)（局部项，已除以 χ_mag）
#   A2 —— 刚度矩阵 A^(2) = (1/4π)[VV − SV − VS + SS]
#         VV 体积散度项（assemble_stiffness!）+ 表面磁荷项 SS−SV−VS（assemble_surface_stiffness!）
#   A  —— 合并 A = A1 + A2（线性系统 A·u = b 的系数矩阵）
#
# 运行方式：
#   julia --project=. src/scripts/build_matrices.jl
#
# 输出（data/ 目录，空格分隔的文本矩阵，可用 readdlm 读回）：
#   data/A1.txt   data/A2.txt   data/A.txt
#
# 说明文件：data/matrices_README.txt

include(joinpath(@__DIR__, "..", "mesh.jl"))
include(joinpath(@__DIR__, "..", "volume.jl"))
include(joinpath(@__DIR__, "..", "mass.jl"))
include(joinpath(@__DIR__, "..", "load.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "surface_stiffness.jl"))
using LinearAlgebra
using DelimitedFiles

function print_status(name, M)
    println("=== ", name, " ===")
    println("  尺寸           = ", size(M))
    println("  sum            = ", sum(M))
    println("  非零元素数     = ", count(!iszero, M))
    println("  对角和 trace   = ", sum(diag(M)))
    println("  对角 min/max   = ", minimum(diag(M)), " / ", maximum(diag(M)))
    println("  行和 max       = ", maximum(abs.(sum(M, dims=2))))
    println("  列和 max       = ", maximum(abs.(sum(M, dims=1))))
    println("  对称性偏差     = ", maximum(abs.(M - M')))
    println("  最小特征值实部 = ", minimum(real.(eigvals(M))))
    println()
end

function main()
    nodes, elements = load_mesh()
    Nn = size(nodes, 1)
    Ne = size(elements, 1)
    ndof = 3 * Nn
    chi_mag = 100.0
    H0 = 1.0
    println("Nn=$Nn  Ne=$Ne  ndof=$ndof  chi_mag=$chi_mag\n")

    # 质量矩阵 A1
    A1 = zeros(ndof, ndof)
    assemble_mass!(A1, elements, nodes, chi_mag)

    # 刚度矩阵 A2 = VV（体积散度项）+ SS − SV − VS（表面磁荷项，完整对称弱形式）
    A2 = zeros(ndof, ndof)
    assemble_stiffness!(A2, elements, nodes; precision=:fast)   # VV 体积项
    assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-3)  # SS − SV − VS 表面项

    # 合并
    A = A1 + A2

    print_status("A1 质量矩阵 A^(1)", A1)
    print_status("A2 刚度矩阵 A^(2)", A2)
    print_status("A = A1 + A2", A)

    # 保存矩阵（空格分隔文本，readdlm 可读回）
    datadir = joinpath(@__DIR__, "..", "..", "data")
    writedlm(joinpath(datadir, "A1.txt"), A1)
    writedlm(joinpath(datadir, "A2.txt"), A2)
    writedlm(joinpath(datadir, "A.txt"), A)
    println("矩阵已保存：data/A1.txt, data/A2.txt, data/A.txt")

    # 说明文件
    open(joinpath(datadir, "matrices_README.txt"), "w") do io
        println(io, "本目录矩阵由 src/scripts/build_matrices.jl 生成。")
        println(io, "运行：julia --project=. src/scripts/build_matrices.jl")
        println(io)
        println(io, "网格：my_box.msh（Nn=", Nn, ", Ne=", Ne, "，ndof=", ndof, "）")
        println(io, "物理参数：chi_mag=", chi_mag, ", H_ext=(H0,0,0), H0=", H0)
        println(io, "积分精度：precision=fast（method=hybrid, rtol=1e-3, atol=1e-6）")
        println(io)
        println(io, "A1.txt  — 质量矩阵 A^(1) = (1/χ_mag)·∫ φ_j·φ_k dΩ（分块对角）")
        println(io, "A2.txt  — 刚度矩阵 A^(2) = (1/4π)[VV − SV − VS + SS]（体积散度 + 表面磁荷，完整对称弱形式）")
        println(io, "A.txt   — 合并 A = A1 + A2，线性系统 A·u = b 的系数矩阵")
        println(io)
        println(io, "读回 Julia：using DelimitedFiles; A = readdlm(\"A.txt\")")
    end
    println("说明已保存：data/matrices_README.txt")
end

main()
