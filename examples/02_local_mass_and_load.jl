# Step 2: understand the mass matrix and load vector from one tetrahedron's local formulas.
#
# 第 2 步：从一个四面体的局部公式理解质量矩阵和载荷向量。

include(joinpath(@__DIR__, "..", "src", "volume.jl"))
include(joinpath(@__DIR__, "..", "src", "mass.jl"))
include(joinpath(@__DIR__, "..", "src", "load.jl"))

nodes = [0.0 0.0 0.0;
         1.0 0.0 0.0;
         0.0 1.0 0.0;
         0.0 0.0 1.0]
elements = reshape([1, 2, 3, 4], 1, 4)

chi_mag = 2.0
H0 = 3.0
ndof = 3 * size(nodes, 1)
A1 = zeros(ndof, ndof)
b = zeros(ndof)

# Julia 约定：函数名末尾的 `!` 表示它会就地修改传入的数组。
assemble_mass!(A1, elements, nodes, chi_mag)
assemble_load!(b, elements, nodes, H0)

# `A1[1:6, 1:6]` 是切片；`1:6` 表示从 1 到 6 的整数范围。
println("A1 左上角 6×6 块：\n", A1[1:6, 1:6])
println("b = ", b)
println("载荷向量总和 = ", sum(b), "，期望 H0*V = ", H0 / 6)
