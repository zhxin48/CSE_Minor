# Topic 3: Aharoni reference quantities and physical checks of the full solution.
#
# 专题 3：Aharoni 参考量和完整解的物理检查。

using DelimitedFiles
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "mesh.jl"))
include(joinpath(@__DIR__, "..", "src", "volume.jl"))
include(joinpath(@__DIR__, "..", "src", "solve.jl"))

chi_mag = 100.0
H0 = 1.0
lx, ly, lz = 1.0, 1.0, 0.1

Nx, Ny, Nz = demag_factors(lx, ly, lz)
uniform_Mx = chi_mag * H0 / (1 + chi_mag * Nx)

println("Aharoni magnetometric demagnetizing factors:")
println("  Nx = ", Nx)
println("  Ny = ", Ny)
println("  Nz = ", Nz)
println("  Nx+Ny+Nz = ", Nx + Ny + Nz, "（应接近 1）")
println("均匀磁化近似 Mx = ", uniform_Mx)

solution_path = joinpath(@__DIR__, "..", "data", "u.txt")
if !isfile(solution_path)
    println("""

尚未找到 data/u.txt，因此这里只完成了解析参考量验证。
先运行下面任一命令生成完整解，再重新运行本专题：

    julia --project=. src/main.jl
    julia --project=. src/scripts/solve_structured.jl
""")
else
    nodes, elements = load_mesh()
    u = vec(readdlm(solution_path, Float64))
    expected_length = 3 * size(nodes, 1)
    length(u) == expected_length || error(
        "u.txt 长度 $(length(u)) 与当前 mesh 的 ndof=$expected_length 不一致")

    total_volume = volume_integral(_ -> 1.0, nodes, elements)
    mean_Mx = volume_integral(i -> Mx(u, i), nodes, elements) / total_volume
    mean_My = volume_integral(i -> My(u, i), nodes, elements) / total_volume
    mean_Mz = volume_integral(i -> Mz(u, i), nodes, elements) / total_volume
    effective_Nx = (chi_mag * H0 / mean_Mx - 1) / chi_mag

    mx = u[1:3:end]
    my = u[2:3:end]
    mz = u[3:3:end]

    println("\n完整离散解检查：")
    println("  mean(Mx,My,Mz) = ", (mean_Mx, mean_My, mean_Mz))
    println("  max|My|/max|Mx| = ", maximum(abs.(my)) / maximum(abs.(mx)))
    println("  max|Mz|/max|Mx| = ", maximum(abs.(mz)) / maximum(abs.(mx)))
    println("  solution-derived effective Nx = ", effective_Nx)
    println("  mean Mx / uniform-Aharoni Mx = ", mean_Mx / uniform_Mx)
    println("""

注意：矩形棱柱的自洽磁化通常不均匀，所以完整解的 mean Mx 不要求等于
均匀磁化假设下的 Aharoni 近似。需要检查的是残差、对称性、网格收敛和
magnetometric factor 的独立积分验证，而不是强迫这两个 Mx 完全相等。
""")
end
