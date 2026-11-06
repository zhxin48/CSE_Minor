# 第 3 步：比较 6D 嵌套积分（方法 A）和 Euler 降维积分（方法 B）。

include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))

reference = ([0.0, 0.0, 0.0],
             [1.0, 0.0, 0.0],
             [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0])

# `x .+ shift` 中的点号表示逐元素相加（broadcasting）。
shift = [5.0, 0.0, 0.0]
separated = tuple((x .+ shift for x in reference)...)

# 分号后的参数是 keyword arguments；调用时可以按名字给出。
method_a = tet_tet_integral_num(
    reference, separated; rtol=1e-3, atol=1e-8, maxevals=100_000)
method_b = tet_tet_integral_facets(
    reference, separated; rtol=1e-3, atol=1e-8, maxevals=100_000)

println("方法 A（6D 嵌套 HCubature） = ", method_a)
println("方法 B（Euler 6D→4D）      = ", method_b)
println("两者绝对差                 = ", abs(method_a - method_b))
