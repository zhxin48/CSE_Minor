# 求解与物理解验。依赖 volume.jl（tet_volume）已先 include。
using LinearAlgebra
using DelimitedFiles

# u 布局：节点 i 对应 Mx = u[3i-2]、My = u[3i-1]、Mz = u[3i]
Mx(u, i) = u[3i - 2]
My(u, i) = u[3i - 1]
Mz(u, i) = u[3i]

"""
    solve_system(A, b) -> u

解线性系统 `A·u = b`（A 对称正定，`A \\ b` 走 LAPACK Cholesky）。
打印条件数与相对残差，返回解向量。
"""
function solve_system(A, b)
    u = A \ b
    rel_res = norm(A * u - b) / norm(b)
    println("求解 A·u = b：")
    println("  cond(A)                 = ", cond(A))
    println("  相对残差 ‖Au-b‖/‖b‖     = ", rel_res)
    return u
end

# P1 基函数上 ∫φ_i dV = V/4，故标量函数 f（节点索引 → 值）的体积分为
#   ∫f dV = Σ_el (V_el/4) Σ_{i∈el} f(i)
function volume_integral(f, nodes, elements)
    total = 0.0
    for alpha in axes(elements, 1)
        inode = elements[alpha, :]
        x = [nodes[inode[k], :] for k in 1:4]
        V = tet_volume(x...)
        total += (V / 4) * sum(f, inode)
    end
    return total
end

"""
    demag_factors(a, b, c) -> (Nx, Ny, Nz)

Aharoni (1998) 矩形棱柱 `a×b×c` 的磁测（magnetometric）退磁因子解析式
（*J. Appl. Phys.* 83, 3432）。返回沿 x/y/z 轴磁化时的退磁因子。

自检恒等式：对任意形状 Nx+Ny+Nz = 1。
"""
function demag_factors(a, b, c)
    Nx = _demag_z(b, c, a)
    Ny = _demag_z(c, a, b)
    Nz = _demag_z(a, b, c)
    return Nx, Ny, Nz
end

# 棱柱 a×b×c 沿 c 轴磁化的退磁因子 Nz（返回 π·Nz，最后除 π）。
function _demag_z(a, b, c)
    s  = sqrt(a^2 + b^2 + c^2)
    t  = 0.0
    t += (b^2 - c^2) / (2*b*c) * log((s - a) / (s + a))
    t += (a^2 - c^2) / (2*a*c) * log((s - b) / (s + b))
    t += b / (2*c) * log((sqrt(a^2 + b^2) + a) / (sqrt(a^2 + b^2) - a))
    t += a / (2*c) * log((sqrt(a^2 + b^2) + b) / (sqrt(a^2 + b^2) - b))
    t += c / (2*a) * log((sqrt(b^2 + c^2) - b) / (sqrt(b^2 + c^2) + b))
    t += c / (2*b) * log((sqrt(a^2 + c^2) - a) / (sqrt(a^2 + c^2) + a))
    t += 2 * atan(a*b / (c*s))
    t += (a^3 + b^3 - 2*c^3) / (3*a*b*c)
    t += (a^2 + b^2 - 2*c^2) / (3*a*b*c) * s
    t += c / (a*b) * (sqrt(a^2 + c^2) + sqrt(b^2 + c^2))
    t -= ((a^2 + b^2)^(3/2) + (b^2 + c^2)^(3/2) + (c^2 + a^2)^(3/2)) / (3*a*b*c)
    return t / π
end

"""
    verify_solution(A, b, u, nodes, elements; H0, chi_mag)

物理解验（假设 H_ext = (H0,0,0)，板关于三个中面镜像对称）：

1. 残差（数值正确性）。
2. 体积分对称性：M_y、M_z 是镜像奇函数 → ∫M_y dV = ∫M_z dV = 0，∫M_x dV > 0。
3. 上界与主导：max|M| ≤ χ_mag·H0；M_x 主导（M_y/M_z 幅度远小于 M_x）。
4. 中心处（质心最近节点）M_y、M_z ≈ 0。
5. 报告平均 M_x，并将 Aharoni 均匀磁化结果作为参考量列出。
"""
function verify_solution(A, b, u, nodes, elements; H0, chi_mag)
    Nn = size(nodes, 1)
    println("\n===== 物理解验 =====")

    # 1. 残差
    println("[1] 相对残差 ‖Au-b‖/‖b‖ = ", norm(A * u - b) / norm(b))

    # 2. 体积分（镜像对称）
    int_Mx = volume_integral(i -> Mx(u, i), nodes, elements)
    int_My = volume_integral(i -> My(u, i), nodes, elements)
    int_Mz = volume_integral(i -> Mz(u, i), nodes, elements)
    println("[2] 体积分：∫Mx dV = ", int_Mx, "，∫My dV = ", int_My, "，∫Mz dV = ", int_Mz)
    println("    对称性要求 ∫My = ∫Mz = 0，∫Mx > 0")

    # 3. 上界与主导
    mx  = [Mx(u, i) for i in 1:Nn]
    my  = [My(u, i) for i in 1:Nn]
    mz  = [Mz(u, i) for i in 1:Nn]
    mn  = [norm([mx[i], my[i], mz[i]]) for i in 1:Nn]
    println("[3] max|M| = ", maximum(mn), " ≤ χ·H0 = ", chi_mag * H0)
    println("    max|Mx| = ", maximum(abs.(mx)),
            "，max|My| = ", maximum(abs.(my)),
            "，max|Mz| = ", maximum(abs.(mz)))
    println("    幅度比 My/Mx = ", maximum(abs.(my)) / maximum(abs.(mx)),
            "，Mz/Mx = ", maximum(abs.(mz)) / maximum(abs.(mx)))

    # 4. 中心处
    c = [(maximum(nodes[:, d]) + minimum(nodes[:, d])) / 2 for d in 1:3]
    ic = argmin([norm(nodes[i, :] - c) for i in 1:Nn])
    println("[4] 中心节点 ", ic, " @ ", nodes[ic, :],
            "：M = (", Mx(u, ic), ", ", My(u, ic), ", ", Mz(u, ic), ")")
    println("    对称性要求中心处 My、Mz ≈ 0")

    # 5. 平均 Mx 与解析退磁因子对比
    V_total = volume_integral(_ -> 1.0, nodes, elements)
    mean_Mx = int_Mx / V_total
    println("[5] 平均 Mx = ", mean_Mx, "（χ·H0 = ", chi_mag * H0, "）")
    Nx_emp = (chi_mag * H0 / mean_Mx - 1) / chi_mag
    println("    经验退磁因子 Nx = ", Nx_emp)

    bbox = [maximum(nodes[:, d]) - minimum(nodes[:, d]) for d in 1:3]
    Nx, Ny, Nz = demag_factors(bbox...)
    if abs(Nx + Ny + Nz - 1) < 1e-3
        Mx_pred = chi_mag * H0 / (1 + chi_mag * Nx)
        println("    Aharoni 解析 Nx = ", Nx, "（Nx+Ny+Nz = ", Nx + Ny + Nz, "）→ 预测 Mx = ", Mx_pred)
        ratio = mean_Mx / Mx_pred
        println("    自洽解平均值/均匀磁化近似 = ", ratio)
        println("    注：矩形棱柱的自洽磁化通常不均匀，因此该比值不要求等于 1")
    else
        @warn "Aharoni 退磁因子自检失败（Nx+Ny+Nz = $(Nx+Ny+Nz) ≠ 1），跳过解析对比"
    end
    return u
end

"""
    save_solution(u, path)

把解向量写为文本（`writedlm`），供可视化等复用。
"""
function save_solution(u, path)
    d = dirname(path)
    isempty(d) || mkpath(d)
    writedlm(path, u)
    println("解已保存：", path)
end
