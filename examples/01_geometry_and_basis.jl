# Step 1: understand a single tetrahedron only — no global mesh or linear system.
#
# 第 1 步：只理解一个四面体，不涉及全局网格或线性系统。

include(joinpath(@__DIR__, "..", "src", "volume.jl"))
include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))

x1 = [0.0, 0.0, 0.0]
x2 = [1.0, 0.0, 0.0]
x3 = [0.0, 1.0, 0.0]
x4 = [0.0, 0.0, 1.0]

# `(x1, x2, x3, x4)` 是 tuple：长度固定，适合表示“恰好四个顶点”。
vertices = (x1, x2, x3, x4)

# `vertices...` 称为 splatting：把 tuple 展开成四个函数参数。
volume = tet_volume(vertices...)
gradients = basis_grads(vertices...)
facets = tet_facets(vertices...)

println("四面体体积 = ", volume, "，期望值 = ", 1 / 6)
println("四个 P1 基函数的梯度（每列对应一个基函数）：\n", gradients)
println("四个面的外法向：")
for (face_number, facet) in enumerate(facets)
    normal = facet[4]
    println("  face ", face_number, ": ", normal)
end
