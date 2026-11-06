# 刚度矩阵核心积分的验证脚本（不集成到 main.jl）
include(joinpath(@__DIR__, "..", "stiffness.jl"))

using LinearAlgebra

# --- 1. integrate_tet 对常数函数应返回参考四面体体积 1/6 ---
x1 = [0.0, 0.0, 0.0]; x2 = [1.0, 0.0, 0.0]
x3 = [0.0, 1.0, 0.0]; x4 = [0.0, 0.0, 1.0]
V = integrate_tet(r -> 1.0, x1, x2, x3, x4)
println("integrate_tet 常数 = ", V, "  (期望 1/6 = ", 1/6, ")")

# --- 2. 两个不相交四面体：参考四面体 + 沿 x 平移 5 ---
ref  = ([0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0])
shift = ([5.0,0.0,0.0], [6.0,0.0,0.0], [5.0,1.0,0.0], [5.0,0.0,1.0])

ref_A = 0.005555531094601571   # notebook 方法 A 参考值
ref_B = 0.005555532350441779   # notebook 方法 B 参考值

println("\n--- 两不相交四面体 I_αβ ---")
IA = tet_tet_integral_num(ref, shift)
IB = tet_tet_integral_facets(ref, shift)
println("方法 A (6D 数值)  = ", IA, "  参考 = ", ref_A, "  差 = ", abs(IA - ref_A))
println("方法 B (4D 降维)  = ", IB, "  参考 = ", ref_B, "  差 = ", abs(IB - ref_B))
println("A 与 B 之差        = ", abs(IA - IB))

# --- 3. 自相互作用（α=β），方法 B 应给出有限值 ---
println("\n--- 自相互作用 I_αα（方法 B）---")
Iself = tet_tet_integral_facets(ref, ref)
println("I_αα (facet) = ", Iself, "  (应为有限正数)")
