# 性能基准：方法 A（6D 数值积分）vs 方法 B（4D 解析降维）
# 扫描不同分离距离 × 不同 rtol，测量耗时与相对误差，用于确定方法 C 的混合策略。
include(joinpath(@__DIR__, "..", "stiffness.jl"))
using Printf

ref_tet = ([0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0])
shift_tet(t) = tuple((ref_tet[i] .+ t for i in 1:4)...)

# 测试用例：沿 x 轴平移的第二个四面体，覆盖自作用→远场
cases = [
    ("自作用 d=0.0",  [0.0, 0.0, 0.0]),
    ("共享顶点 d=1.0",[1.0, 0.0, 0.0]),
    ("近场 d=1.2",    [1.2, 0.0, 0.0]),
    ("近场 d=1.5",    [1.5, 0.0, 0.0]),
    ("中距 d=2.0",    [2.0, 0.0, 0.0]),
    ("中距 d=3.0",    [3.0, 0.0, 0.0]),
    ("远场 d=5.0",    [5.0, 0.0, 0.0]),
]

rtols = [1e-2, 1e-3, 1e-4, 1e-6]
ref_rtol = 1e-7

# 热身（触发 JIT 编译）
println("热身编译中...")
tet_tet_integral_facets(ref_tet, shift_tet([5.0,0,0]); rtol=1e-3)
tet_tet_integral_facets(ref_tet, ref_tet; rtol=1e-3)
tet_tet_integral_num(ref_tet, shift_tet([5.0,0,0]); rtol=1e-3)

function timed(f; n=3)
    ts = Float64[]
    for _ in 1:n
        GC.gc()
        push!(ts, @elapsed f())
    end
    minimum(ts)
end

# 参考值：方法 B 高精度（对所有用例都可靠，含自作用）
println("计算参考值（方法 B, rtol=$ref_rtol）...")
Irefs = Dict{String,Float64}()
for (name, t) in cases
    B = shift_tet(t)
    Irefs[name] = tet_tet_integral_facets(ref_tet, B; rtol=ref_rtol, atol=1e-9)
    println("  $name -> ", Irefs[name])
end

# 测量
results = []
for (name, t) in cases
    B = shift_tet(t)
    Iref = Irefs[name]
    d = norm(t)
    for rtol in rtols
        # 方法 B
        vB = tet_tet_integral_facets(ref_tet, B; rtol=rtol, atol=1e-6)
        tB = timed(() -> tet_tet_integral_facets(ref_tet, B; rtol=rtol, atol=1e-6))
        push!(results, (name, d, "B", rtol, vB, abs(vB-Iref)/abs(Iref), tB*1000))

        # 方法 A（自作用奇异，跳过）
        if d > 0.0
            vA = tet_tet_integral_num(ref_tet, B; rtol=rtol, atol=1e-6, maxevals=10^6)
            tA = timed(() -> tet_tet_integral_num(ref_tet, B; rtol=rtol, atol=1e-6, maxevals=10^6))
            push!(results, (name, d, "A", rtol, vA, abs(vA-Iref)/abs(Iref), tA*1000))
        end
        println("  完成 $name rtol=$rtol")
        flush(stdout)
    end
end

# 输出 CSV
out = joinpath(@__DIR__, "..", "..", "data", "benchmark.csv")
open(out, "w") do io
    println(io, "case,distance,method,rtol,value,rel_err,time_ms")
    for r in results
        println(io, join(r, ","))
    end
end

# 打印表格
println("\ncase               d     method  rtol    rel_err      time_ms")
for r in results
    name, d, m, rtol, v, err, tms = r
    @printf("%-16s %5.1f   %s    %.0e  %.3e  %8.2f\n", name, d, m, rtol, err, tms)
end
println("\n结果已写入 ", out)
