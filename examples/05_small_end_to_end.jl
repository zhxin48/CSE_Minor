# 第 5 步：用单个四面体走完整流程，避免第一次复现就等待完整板网格。

include(joinpath(@__DIR__, "..", "src", "volume.jl"))
include(joinpath(@__DIR__, "..", "src", "mass.jl"))
include(joinpath(@__DIR__, "..", "src", "load.jl"))
include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "src", "surface_stiffness.jl"))
include(joinpath(@__DIR__, "..", "src", "solve.jl"))

nodes = [0.0 0.0 0.0;
         1.0 0.0 0.0;
         0.0 1.0 0.0;
         0.0 0.0 1.0]
elements = reshape([1, 2, 3, 4], 1, 4)
chi_mag = 2.0
H0 = 1.0
ndof = 3 * size(nodes, 1)

A1 = zeros(ndof, ndof)
A2 = zeros(ndof, ndof)
b = zeros(ndof)

assemble_mass!(A1, elements, nodes, chi_mag)
assemble_stiffness!(A2, elements, nodes;
                    method=:facet, rtol=1e-2, atol=1e-8, maxevals=100_000)
assemble_surface_stiffness!(A2, elements, nodes;
                            rtol=1e-2, atol=1e-8, maxevals=100_000)
assemble_load!(b, elements, nodes, H0)

A = A1 + A2
u = solve_system(A, b)

println("矩阵尺寸 = ", size(A))
println("对称性误差 = ", maximum(abs.(A - A')))
println("相对残差 = ", norm(A * u - b) / norm(b))
println("每个节点的 (Mx, My, Mz)：")
for node in 1:size(nodes, 1)
    println("  node ", node, ": ", u[3node-2:3node])
end
