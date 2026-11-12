# Volume stiffness assembly (the VV volume-divergence term) plus the shared
# quadrature and basis helpers used by surface_stiffness.jl.
using HCubature
using LinearAlgebra

# Simple progress bar (no third-party dependency): refresh every ~1%, newline at i == n.
# 简单进度条（无第三方依赖）：每 ~1% 刷新一次，i==n 时换行。
function _progress(i, n, label; width=40)
    if i == n || (n > 0 && i % max(1, n ÷ 100) == 0)
        frac = i / n
        filled = round(Int, frac * width)
        bar = '#'^filled * '-'^(width - filled)
        print("\r[", label, "] [", bar, "] ", round(100 * frac, digits=1), "%  (", i, "/", n, ")")
        flush(stdout)
        i == n && println()
    end
end

# Gradients of the four P1 basis functions of a tetrahedron (3×4 matrix).
"""
    basis_grads(x1, x2, x3, x4)

返回四面体 4 个 P1 基函数的梯度，3×4 矩阵，第 i 列是 ∇φ_i。
基函数 φ_i 在顶点 i 为 1、其余顶点为 0，在每个四面体内线性。
"""
function basis_grads(x1, x2, x3, x4)
    J = hcat(x2 - x1, x3 - x1, x4 - x1)   # 3×3 仿射映射：x = x1 + J·ξ
    Jinv = inv(J)
    # 参考四面体基函数 φ̂ = (1-ξ-η-ζ, ξ, η, ζ) 对 ξ 的梯度（3×4）
    grad_ref = [-1.0 1.0 0.0 0.0;
                -1.0 0.0 1.0 0.0;
                -1.0 0.0 0.0 1.0]
    return Jinv' * grad_ref               # ∇_x φ = J^{-T} ∇_ξ φ̂
end

# Numerically integrate f over a tetrahedron via hcubature (affine map to the reference tetrahedron).
"""
    integrate_tet(f, x1, x2, x3, x4; ...)

在四面体上做 hcubature 数值积分。
一般四面体 → 参考四面体（仿射，|det J|），再映射到单位立方体 [0,1]^3。
"""
function integrate_tet(f, x1, x2, x3, x4;
                       rtol::Real=1e-5, atol::Real=1e-8,
                       maxevals::Integer=10^7)
    J = hcat(x2 - x1, x3 - x1, x4 - x1)
    detJ = abs(det(J))
    integrand = function (uvw)
        u, v, w = uvw[1], uvw[2], uvw[3]
        jac = (1 - u)^2 * (1 - v)                       # 立方体→参考四面体的 Jacobian
        xi = [u, v * (1 - u), w * (1 - u) * (1 - v)]
        r = x1 + J * xi
        return f(r) * detJ * jac
    end
    # Tuple 端点把积分维数编码进类型；Vector 端点的长度只能在运行时获知，
    # 会使 HCubature 的返回类型无法被 Julia 编译器推断。
    val, _ = hcubature(integrand, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0);
                       rtol=rtol, atol=atol, maxevals=maxevals)
    return val
end

# Numerically integrate f over a triangular facet via hcubature.
"""
    integrate_tri(f, p1, p2, p3; ...)

在三角 facet 上做 hcubature 数值积分（方法 B 的 4D 降维用）。
映射 r(u,v) = p1 + u(p2-p1) + uv(p3-p2)，面积元 2·area·u。
"""
function integrate_tri(f, p1, p2, p3;
                       rtol::Real=1e-5, atol::Real=1e-8,
                       maxevals::Integer=10^7)
    v1 = p2 - p1
    v2 = p3 - p1
    twice_area = norm(cross(v1, v2))
    integrand = function (uv)
        u, v = uv[1], uv[2]
        r = p1 + u * v1 + u * v * (p3 - p2)
        return f(r) * twice_area * u
    end
    val, _ = hcubature(integrand, (0.0, 0.0), (1.0, 1.0);
                       rtol=rtol, atol=atol, maxevals=maxevals)
    return val
end

# Return the four triangular faces of a tetrahedron, each (p1, p2, p3, unit outward normal).
"""
    tet_facets(x1, x2, x3, x4)

返回四面体的 4 个三角面，每个面为 (p1, p2, p3, 单位外法向)。
外法向用"第四顶点在面内侧"判定，与顶点顺序无关。
"""
# Orient facet (p1, p2, p3) so its unit normal points away from `opposite_vertex`.
function _oriented_facet(p1, p2, p3, opposite_vertex)
    n = cross(p2 - p1, p3 - p1)
    n = n / norm(n)
    if dot(n, p1 - opposite_vertex) < 0
        n = -n
    end
    return (p1, p2, p3, n)
end

# Return the four oriented triangular faces of the tetrahedron.
function tet_facets(x1, x2, x3, x4)
    # 显式列出四个面，避免运行期索引使 tuple 的元素类型退化为 Any。
    return (_oriented_facet(x2, x3, x4, x1),
            _oriented_facet(x1, x4, x3, x2),
            _oriented_facet(x1, x2, x4, x3),
            _oriented_facet(x1, x3, x2, x4))
end

# Method A: 6D numerical integration of the Coulomb double volume integral.
"""
    tet_tet_integral_num(alpha_verts, beta_verts; ...)

方法 A：嵌套 hcubature 的 6D 数值积分。
I_{αβ} = ∫_{P_α} ∫_{P_β} 1/||r-r'|| dΩ' dΩ。
仅适用于不相交（远场）四面体——α=β 时核奇异，hcubature 不收敛。
"""
function tet_tet_integral_num(alpha_verts, beta_verts;
                              rtol::Real=1e-5, atol::Real=1e-8,
                              maxevals::Integer=10^7)
    inner = function (rd)                        # 固定 rd ∈ P_α，对 P_β 积分
        value = integrate_tet(rs -> 1 / norm(rd - rs), beta_verts...;
                              rtol=rtol, atol=atol, maxevals=maxevals)
        return value::Float64
    end
    return integrate_tet(inner, alpha_verts...;
                         rtol=rtol, atol=atol, maxevals=maxevals)::Float64
end

# Method B: analytic reduction (Euler homogeneous-function theorem) to 16 4D facet-facet integrals.
"""
    tet_tet_integral_facets(alpha_verts, beta_verts; ...)

方法 B：解析降维（Euler 齐次函数定理）。
6D 体积-体积积分 → 4×4 = 16 个 4D facet-facet 积分：

    I_{αβ} = -1/6 Σ_{F_a} Σ_{F_b} ∫∫ [n_a·(r-r')][n_b·(r-r')]/||r-r'|| dS' dS

被积函数在 r→r' 时是 O(δ)，核不自奇，因此可处理近场与自相互作用 α=β。
"""
function tet_tet_integral_facets(alpha_verts, beta_verts;
                                 rtol::Real=1e-5, atol::Real=1e-8,
                                 maxevals::Integer=10^7)
    facA = tet_facets(alpha_verts...)
    facB = tet_facets(beta_verts...)
    total = 0.0
    for (p1a, p2a, p3a, na) in facA
        for (p1b, p2b, p3b, nb) in facB
            inner = function (rd)                # 固定 rd ∈ F_a，对 F_b 积分
                value = integrate_tri(
                    rs -> dot(na, rd - rs) * dot(nb, rd - rs) / norm(rd - rs),
                    p1b, p2b, p3b;
                    rtol=rtol, atol=atol, maxevals=maxevals)
                return value::Float64
            end
            total += integrate_tri(inner, p1a, p2a, p3a;
                                   rtol=rtol, atol=atol, maxevals=maxevals)
        end
    end
    return -total / 6
end

# Map the `precision` knob to a (method, rtol, atol) triple.
"""
    _stiffness_preset(precision)

「精度 + 近场分派」的单一旋钮 → 返回 (method, rtol, atol)。

precision=:fast     档 1（默认）：自作用 α=β → B，其余 → A；rtol=1e-3（矩阵元素 ~1e-4，快）
precision=:accurate 档 2：自作用 + 共享节点近邻 → B，其余 → A；rtol=1e-4（~1e-5，慢，见 PERFORMANCE_ANALYSIS.md）
"""
function _stiffness_preset(precision)
    precision == :fast     && return (:hybrid,  1e-3, 1e-6)
    precision == :accurate && return (:hybrid2, 1e-4, 1e-6)
    error("未知 precision: $precision（可选 :fast / :accurate）")
end

# Return whether tetrahedra alpha and beta share at least one node.
"""
    share_node(elements, alpha, beta)

四面体 α、β 是否共享至少一个节点（近邻判定，O(1)）。
"""
share_node(elements, alpha, beta) = !isempty(intersect(elements[alpha, :], elements[beta, :]))

# Assemble the global volume stiffness matrix A^(2) (3Nn × 3Nn) into A2 in place.
"""
    assemble_stiffness!(A2, elements, nodes; precision=:fast, ...)

组装全局刚度矩阵 A^(2)（3Nn × 3Nn）。

主接口是 `precision`（单一参数，同时决定积分精度与近场分派）：
    precision=:fast       自作用→B，其余→A，rtol=1e-3
    precision=:accurate   自作用 + 共享节点近邻→B，其余→A，rtol=1e-4

底层 method 仍可单独指定（会覆盖 precision 的预设）：
    method=:facet    全部 B（解析降维 4D，支持近场与自作用）
    method=:num      全部 A（6D 数值，仅远场，α=β 跳过并告警）
    method=:hybrid   自作用→B，其余→A
    method=:hybrid2  自作用 + 共享节点近邻→B，其余→A

局部贡献：Kloc[3(i-1)+a, 3(j-1)+b] = (1/4π) ∇φ_i[a] · ∇φ_j[b] · I_{αβ}
"""
function assemble_stiffness!(A2, elements, nodes; precision=:fast, method=nothing, rtol=nothing, atol=nothing, kws...)
    m0, r0, a0 = _stiffness_preset(precision)
    method = method === nothing ? m0 : method
    rtol   = rtol   === nothing ? r0 : rtol
    atol   = atol   === nothing ? a0 : atol

    Ne = size(elements, 1)
    # I_{alpha,beta} = I_{beta,alpha}. 只计算无序单元对一次，再把转置贡献
    # 写入矩阵另一半。这样既少做一次自适应积分，也保证离散双线性型严格对称，
    # 不会让两次独立积分的误差表现为 A2 != A2'.
    total = Ne * (Ne + 1) ÷ 2
    iter = 0
    for alpha in 1:Ne, beta in alpha:Ne
        iter += 1
        _progress(iter, total, "A^(2) 体积项 VV")
        # 取 α、β 的 4 个顶点坐标
        inode_a = elements[alpha, :]
        inode_b = elements[beta, :]
        xA = [nodes[inode_a[i], :] for i in 1:4]
        xB = [nodes[inode_b[i], :] for i in 1:4]

        # 分派：该单元对用方法 B（解析降维）还是方法 A（6D 数值）
        if method == :facet
            use_facet = true
        elseif method == :num
            if alpha == beta
                @warn "method=:num 无法处理自相互作用 (alpha=$alpha)，跳过"
                continue
            end
            use_facet = false
        elseif method == :hybrid
            use_facet = (alpha == beta)
        elseif method == :hybrid2
            use_facet = (alpha == beta) || share_node(elements, alpha, beta)
        else
            error("未知 method: $method（可选 :facet / :num / :hybrid / :hybrid2）")
        end

        Iab = use_facet ? tet_tet_integral_facets(xA, xB; rtol=rtol, atol=atol, kws...) :
                          tet_tet_integral_num(xA, xB; rtol=rtol, atol=atol, kws...)

        gradsA = basis_grads(xA...)
        gradsB = basis_grads(xB...)

        # 局部 12×12 矩阵
        Kloc = zeros(12, 12)
        for i in 1:4, j in 1:4, a in 1:3, b in 1:3
            gi = 3 * (i - 1) + a
            gj = 3 * (j - 1) + b
            Kloc[gi, gj] = (1 / (4π)) * gradsA[a, i] * gradsB[b, j] * Iab
        end

        # 局部自由度 → 全局自由度映射（节点 i 对应 3i-2/3i-1/3i 为 x/y/z）
        gdof_a = [3 * inode_a[i] - 3 + a for i in 1:4 for a in 1:3]
        gdof_b = [3 * inode_b[j] - 3 + b for j in 1:4 for b in 1:3]
        for ii in 1:12, jj in 1:12
            A2[gdof_a[ii], gdof_b[jj]] += Kloc[ii, jj]
            if alpha != beta
                A2[gdof_b[jj], gdof_a[ii]] += Kloc[ii, jj]
            end
        end
    end
    return A2
end
