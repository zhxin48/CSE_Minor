# Surface demagnetisation terms (Galerkin weak form of the surface magnetic
# charge). Implements SS and the SV/VS cross terms that complete the volume
# divergence term VV from stiffness.jl into A^(2) = (1/4π)[VV − SV − VS + SS].
#
# 表面退磁项（表面磁荷的 Galerkin 弱形式）。
#
# 完整对称弱形式的退磁刚度为
#     A^(2) = (1/4π)[ VV − SV − VS + SS ]
# 其中 VV 是已有的体积-体积散度项（stiffness.jl）。本文件实现表面磁荷项 SS
# 与体积-表面交叉项 SV/VS（VS = SVᵀ）。均匀磁化时 VV=VS=0，退磁由 −SV + SS 给出。
#
# 依赖 stiffness.jl 的 tet_facets / basis_grads；入口脚本应先 include stiffness.jl。
#
# SS 自面积分用 **Duffy 变换**（self_facet_matrix）消除 r=r' 的 2D 奇异对角线；
# 非自面用 4D hcubature（facet_facet_matrix）。SV 用已验证的体积→面 Euler 降维（因子 1/2）。
# 已验证：单个参考四面体 Nxx=1/3（解析）、立方体 trace=1.0，均精确吻合。
using HCubature
using LinearAlgebra

# Find the triangular faces on the mesh boundary (owned by a single tetrahedron).
"""
    boundary_faces(elements, nodes) -> [(i, j, k, n), ...]

找出网格的边界面（只属于一个四面体的三角面）。返回每个面的
三个全局节点索引 (i,j,k) 与单位外法向 n。
"""
function boundary_faces(elements, nodes)
    count = Dict{NTuple{3,Int},Int}()
    owner = Dict{NTuple{3,Int},Int}()
    orig  = Dict{NTuple{3,Int},NTuple{3,Int}}()
    for alpha in axes(elements, 1)
        v = elements[alpha, :]
        for (a, b, c) in ((2, 3, 4), (1, 4, 3), (1, 2, 4), (1, 3, 2))
            key = Tuple(sort([v[a], v[b], v[c]]))
            count[key] = get(count, key, 0) + 1
            owner[key] = alpha
            orig[key]  = (v[a], v[b], v[c])
        end
    end
    faces = Vector{Tuple{Int,Int,Int,Vector{Float64}}}()
    for (key, c) in count
        c == 1 || continue
        i, j, k = orig[key]
        alpha = owner[key]
        x1, x2, x3 = nodes[i, :], nodes[j, :], nodes[k, :]
        d = first(setdiff(elements[alpha, :], [i, j, k]))
        x4 = nodes[d, :]
        n = cross(x2 - x1, x3 - x1)
        n = n / norm(n)
        if dot(n, x1 - x4) < 0   # 外法向：指向远离内部顶点 x4 一侧
            n = -n
        end
        push!(faces, (i, j, k, n))
    end
    return faces
end

# Single-layer potential matrix between two triangular facets (4D hcubature).
"""
    facet_facet_matrix(fA, fB; rtol, atol, maxevals) -> 3×3 matrix

两个三角面 fA=(a1,a2,a3)、fB=(b1,b2,b3) 之间的单层位势矩阵
    S[m,n] = ∫_{fA} ∫_{fB} φ_m(r) φ_n(r') / |r − r'| dS' dS
其中 φ_m 是面上第 m 个 P1 基（重心坐标）。4D hcubature 直接求积；
1/|r−r'| 奇异性可积（自面时在 2D 对角线上，雅可比抵消之）。
"""
function facet_facet_matrix(fA, fB;
                            rtol::Real=1e-4, atol::Real=1e-8,
                            maxevals::Integer=10^7)
    a1, a2, a3 = fA
    b1, b2, b3 = fB
    ja = norm(cross(a2 - a1, a3 - a1))   # 2·area(A)
    jb = norm(cross(b2 - b1, b3 - b1))   # 2·area(B)
    e23a = a3 - a2
    e23b = b3 - b2
    integrand = function (uv)
        uA, vA, uB, vB = uv[1], uv[2], uv[3], uv[4]
        rA = a1 + uA * (a2 - a1) + uA * vA * e23a
        rB = b1 + uB * (b2 - b1) + uB * vB * e23b
        phiA = [1 - uA, uA * (1 - vA), uA * vA]
        phiB = [1 - uB, uB * (1 - vB), uB * vB]
        return (ja * jb * uA * uB / norm(rA - rB)) .* (phiA * phiB')
    end
    S, _ = hcubature(integrand, (0.0, 0.0, 0.0, 0.0),
                     (1.0, 1.0, 1.0, 1.0);
                     rtol=rtol, atol=atol, maxevals=maxevals)
    return S::Matrix{Float64}
end

# Signed barycentric coordinates of point r in the triangle (a1, a2, a3).
"""
    tri_barycentric(a1, a2, a3, r) -> 3-vector

点 r 在三角面 (a1,a2,a3) 内的重心坐标（带符号，满足 Σλ = 1）。
"""
function tri_barycentric(a1, a2, a3, r)
    n = cross(a2 - a1, a3 - a1)
    n2 = dot(n, n)
    return [dot(cross(a2 - r, a3 - r), n),
            dot(cross(a3 - r, a1 - r), n),
            dot(cross(a1 - r, a2 - r), n)] ./ n2
end

# Self single-layer potential matrix using a Duffy transform to remove the singular diagonal.
"""
    self_facet_matrix(fF; kws...) -> 3×3 matrix

自面单层位势 S[m,n] = ∫_F ∫_F φ_m(r) φ_n(r') / |r−r'| dS' dS，用 **Duffy 变换**
消除 r=r' 的奇异对角线：外层 r 用常规三角映射，内层把 F 沿 r 切成 3 个子三角，
每个子三角内用 Duffy 映射 r' = r + u(va−r) + uv(vb−va)，使 dS'/|r−r'| 变为
非奇异 (jac/denom)。与边×边 Euler 归约交叉验证一致（参考三角 sum(S)=1.00307）。
"""
function self_facet_matrix(fF;
                           rtol::Real=1e-4, atol::Real=1e-8,
                           maxevals::Integer=10^7)
    a1, a2, a3 = fF
    inner = function (r)
        total = zeros(3)
        for (va, vb) in ((a1, a2), (a2, a3), (a3, a1))
            jac = norm(cross(va - r, vb - r))
            var = va - r
            vbva = vb - va
            integ = function (uv)
                u, v = uv[1], uv[2]
                rp = r + u * var + u * v * vbva
                denom = norm(var + v * vbva)
                phi = tri_barycentric(a1, a2, a3, rp)
                return (jac / denom) .* phi
            end
            total .+= first(hcubature(integ, (0.0, 0.0), (1.0, 1.0); rtol=rtol, atol=atol, maxevals=maxevals))
        end
        return total::Vector{Float64}
    end
    ja = norm(cross(a2 - a1, a3 - a1))
    e23 = a3 - a2
    integ = function (uv)
        u, v = uv[1], uv[2]
        r = a1 + u * (a2 - a1) + u * v * e23
        phi = [1 - u, u * (1 - v), u * v]
        return (ja * u) .* (phi * transpose(inner(r)))
    end
    S, _ = hcubature(integ, (0.0, 0.0), (1.0, 1.0); rtol=rtol, atol=atol, maxevals=maxevals)
    return S::Matrix{Float64}
end

# Surface-volume cross-term vector W_j = ∫∫ φ_j (r−r')·nP / |r−r'| dS dS'.
"""
    facet_facet_sv_vector(fF, fP, nP; kws...) -> 3-vector

计算 W_j = ∫_{fF} ∫_{fP} φ_j(r') (r − r')·nP / |r − r'| dS_r dS_{r'}，
其中 φ_j 是面 fF 上第 j 个 P1 基（r' 一侧），nP 是面 fP 的外法向。
核 (r−r')·nP/|r−r'| 有界（余弦），非奇异，hcubature 直接求积。
"""
function facet_facet_sv_vector(fF, fP, nP;
                               rtol::Real=1e-4, atol::Real=1e-8,
                               maxevals::Integer=10^7)
    a1, a2, a3 = fF
    b1, b2, b3 = fP
    ja = norm(cross(a2 - a1, a3 - a1))
    jb = norm(cross(b2 - b1, b3 - b1))
    e23a = a3 - a2
    e23b = b3 - b2
    integrand = function (uv)
        uA, vA, uB, vB = uv[1], uv[2], uv[3], uv[4]
        rA = a1 + uA * (a2 - a1) + uA * vA * e23a   # r' ∈ fF
        rB = b1 + uB * (b2 - b1) + uB * vB * e23b   # r  ∈ fP
        phi = [1 - uA, uA * (1 - vA), uA * vA]
        R = rB - rA
        w = (ja * jb * uA * uB) * dot(R, nP) / norm(R)
        return w .* phi
    end
    W, _ = hcubature(integrand, (0.0, 0.0, 0.0, 0.0),
                     (1.0, 1.0, 1.0, 1.0);
                     rtol=rtol, atol=atol, maxevals=maxevals)
    return W::Vector{Float64}
end

# Accumulate the surface demagnetisation terms into A2 in place:
# A2 += (1/4π)SS − (1/4π)(SV + SVᵀ).
"""
    assemble_surface_stiffness!(A2, elements, nodes; rtol, atol, maxevals)

把表面退磁项累加进 A2（就地）：
    A2 += (1/4π) SS − (1/4π)(SV + SVᵀ)
即补上体积散度项 VV（stiffness.jl）所缺的表面磁荷部分。
"""
function assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-4, atol=1e-8, maxevals=10^7)
    bfaces = boundary_faces(elements, nodes)
    Nf = length(bfaces)
    Ne = size(elements, 1)
    inv4pi = 1 / (4π)

    # 预计算：每个四面体的 4 个 facet（含外法向）与基函数梯度
    tet_fac = [tet_facets([nodes[elements[alpha, q], :] for q in 1:4]...) for alpha in 1:Ne]
    tet_gr  = [basis_grads([nodes[elements[alpha, q], :] for q in 1:4]...) for alpha in 1:Ne]
    # 预计算：每个边界面的坐标 / 节点 / 法向
    face_coords = [(nodes[bfaces[a][1], :], nodes[bfaces[a][2], :], nodes[bfaces[a][3], :]) for a in 1:Nf]
    face_nodes  = [(bfaces[a][1], bfaces[a][2], bfaces[a][3]) for a in 1:Nf]
    face_n      = [bfaces[a][4] for a in 1:Nf]

    # --- SS：面-面单层位势 ---
    # S_BA = S_AB'. 只积分无序面面对，然后显式装配转置块，避免两次
    # 自适应积分产生略有不同的结果而破坏全局矩阵的理论对称性。
    total_ss = Nf * (Nf + 1) ÷ 2
    iter = 0
    for a in 1:Nf, b in a:Nf
        iter += 1
        _progress(iter, total_ss, "表面项 SS")
        S = if a == b
            self_facet_matrix(face_coords[a]; rtol=rtol, atol=atol, maxevals=maxevals)
        else
            facet_facet_matrix(face_coords[a], face_coords[b]; rtol=rtol, atol=atol, maxevals=maxevals)
        end
        a == b && (S = (S + S') / 2)
        for m in 1:3, n in 1:3
            gA = face_nodes[a][m]
            gB = face_nodes[b][n]
            for cA in 1:3, cB in 1:3
                row = 3*(gA-1)+cA
                col = 3*(gB-1)+cB
                val = inv4pi * face_n[a][cA] * face_n[b][cB] * S[m, n]
                A2[row, col] += val
                if a != b
                    A2[col, row] += val
                end
            end
        end
    end

    # --- SV + VS（= SVᵀ）：面电荷 × 体散度交叉项，负号 ---
    total_sv = Nf * Ne
    iter = 0
    for a in 1:Nf, alpha in 1:Ne
        iter += 1
        _progress(iter, total_sv, "表面项 SV/VS")
        # V_j = ∫_F ∫_P φ_j / R = (1/2) Σ_{f⊂P} ∫_F ∫_f φ_j (r−r')·n_f / |r−r'|
        V = zeros(3)
        for (p1, p2, p3, nf) in tet_fac[alpha]
            V .+= facet_facet_sv_vector(face_coords[a], (p1, p2, p3), nf; rtol=rtol, atol=atol, maxevals=maxevals)
        end
        V .*= 0.5

        gradsA = tet_gr[alpha]            # gradsA[aa, l] = ∂_aa φ_{node l}
        inode = elements[alpha, :]
        for m in 1:3, l in 1:4, b in 1:3, aa in 1:3
            val = inv4pi * gradsA[aa, l] * face_n[a][b] * V[m]
            k = inode[l]                  # 测试节点
            j = face_nodes[a][m]          # 试节点（面上）
            A2[3*(k-1)+aa, 3*(j-1)+b] -= val   # SV
            A2[3*(j-1)+b, 3*(k-1)+aa] -= val   # VS = SVᵀ
        end
    end
    return A2
end
