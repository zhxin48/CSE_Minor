# Mass-matrix (A^(1)) assembly for the magnetisation weak form.
using LinearAlgebra

# Local scalar mass matrix of a 4-node tetrahedron.
"""
    local_mass(V)

4 节点四面体的局部标量质量矩阵。
"""
function local_mass(V)
    Mloc = fill(V/20, 4, 4)
    for i in 1:4
        Mloc[i, i] = V/10
    end
    return Mloc
end

# Assemble the global mass matrix A1 (3Nn × 3Nn), scaled by 1/chi_mag.
"""
    assemble_mass!(A1, elements, nodes, chi_mag)

组装全局质量矩阵 A1（3Nn × 3Nn），并包含系数 `1/chi_mag`。

函数名末尾的 `!` 是 Julia 约定：函数会就地修改第一个参数 `A1`。
"""
function assemble_mass!(A1, elements, nodes, chi_mag)
    Ne = size(elements, 1)
    for alpha in 1:Ne
        inode = elements[alpha, :]
        x1 = nodes[inode[1], :]
        x2 = nodes[inode[2], :]
        x3 = nodes[inode[3], :]
        x4 = nodes[inode[4], :]
        V = tet_volume(x1, x2, x3, x4)
        Mloc = local_mass(V)

        # 局部 12 个自由度对应的全局自由度编号
        gdof = zeros(Int, 12)
        for i in 1:4
            gdof[3*(i-1) + 1] = 3*inode[i] - 2   # x
            gdof[3*(i-1) + 2] = 3*inode[i] - 1   # y
            gdof[3*(i-1) + 3] = 3*inode[i]       # z
        end

        # 累加：分量 a 和 b 相同才有贡献
        for a in 1:3
            for b in 1:3
                if a != b
                    continue
                end
                for i in 1:4
                    for j in 1:4
                        gi = gdof[3*(i-1) + a]
                        gj = gdof[3*(j-1) + b]
                        A1[gi, gj] += Mloc[i, j] / chi_mag
                    end
                end
            end
        end
    end
end
