# 专题 1：多单元 scatter，以及自作用/近场/远场的 A-B 混合分派。

include(joinpath(@__DIR__, "..", "src", "volume.jl"))
include(joinpath(@__DIR__, "..", "src", "mass.jl"))
include(joinpath(@__DIR__, "..", "src", "load.jl"))
include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))

# 前两个四面体共享三角面 (2,3,4)。后四个节点组成一个远处四面体。
nodes = [0.0 0.0 0.0;   # 1
         1.0 0.0 0.0;   # 2
         0.0 1.0 0.0;   # 3
         0.0 0.0 1.0;   # 4
         1.0 1.0 1.0;   # 5
         5.0 0.0 0.0;   # 6
         6.0 0.0 0.0;   # 7
         5.0 1.0 0.0;   # 8
         5.0 0.0 1.0]   # 9

adjacent_elements = [1 2 3 4;
                     5 2 4 3]

ndof = 3 * size(nodes, 1)
A1 = zeros(ndof, ndof)
b = zeros(ndof)
assemble_mass!(A1, adjacent_elements, nodes, 2.0)
assemble_load!(b, adjacent_elements, nodes, 1.0)

println("共享面节点 2、3、4 会同时收到两个单元的局部贡献。")
for node in (1, 2, 3, 4, 5)
    x_dof = 3 * node - 2
    println("  node ", node,
            "  A1[x,x]=", A1[x_dof, x_dof],
            "  b[x]=", b[x_dof])
end

# 与 assemble_stiffness! 的 hybrid2 判据一致：
# 自作用和共享节点近场用 B；完全分离的远场用 A。
function selected_method(elements, alpha, beta)
    if alpha == beta
        return :B_self
    elseif share_node(elements, alpha, beta)
        return :B_near
    else
        return :A_far
    end
end

all_elements = [1 2 3 4;
                5 2 4 3;
                6 7 8 9]
println("\n混合分派：")
for (alpha, beta) in ((1, 1), (1, 2), (1, 3))
    println("  pair (", alpha, ",", beta, ") → ",
            selected_method(all_elements, alpha, beta))
end

tet1 = tuple((nodes[i, :] for i in all_elements[1, :])...)
tet2 = tuple((nodes[i, :] for i in all_elements[2, :])...)
tet3 = tuple((nodes[i, :] for i in all_elements[3, :])...)

# 只计算三种代表关系，而不是组装全部 3×3 个单元对。
self_value = tet_tet_integral_facets(
    tet1, tet1; rtol=1e-2, atol=1e-8, maxevals=100_000)
near_value = tet_tet_integral_facets(
    tet1, tet2; rtol=1e-2, atol=1e-8, maxevals=100_000)
far_value = tet_tet_integral_num(
    tet1, tet3; rtol=1e-2, atol=1e-8, maxevals=100_000)

println("\n代表积分：")
println("  self，方法 B = ", self_value)
println("  near，方法 B = ", near_value)
println("  far， 方法 A = ", far_value)
