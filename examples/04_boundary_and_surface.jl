# Step 4: identify boundary faces and compute the Duffy self-interaction matrix of one face.
#
# 第 4 步：识别边界面，并计算一个面的 Duffy 自相互作用矩阵。

include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "src", "surface_stiffness.jl"))

nodes = [0.0 0.0 0.0;
         1.0 0.0 0.0;
         0.0 1.0 0.0;
         0.0 0.0 1.0]
elements = reshape([1, 2, 3, 4], 1, 4)

faces = boundary_faces(elements, nodes)
println("单个四面体共有 ", length(faces), " 个边界面。")
for (number, (i, j, k, normal)) in enumerate(faces)
    println("  face ", number, " nodes=", (i, j, k), " normal=", normal)
end

i, j, k, _ = first(faces)
face_coordinates = (nodes[i, :], nodes[j, :], nodes[k, :])
self_matrix = self_facet_matrix(
    face_coordinates; rtol=1e-3, atol=1e-8, maxevals=100_000)

println("第一个面的 Duffy 自相互作用矩阵：\n", self_matrix)
println("对称性误差 = ", maximum(abs.(self_matrix - self_matrix')))
