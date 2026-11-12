# Load-vector assembly for the uniform external magnetic field.
# Assemble the global load vector b (length 3*Nn) for H_ext = (H0, 0, 0).
"""
    assemble_load!(b, elements, nodes, H0)

组装全局载荷向量 b（长度 3Nn）。
假设外部磁场 H_ext = (H0, 0, 0)。
"""
function assemble_load!(b, elements, nodes, H0)
    Ne = size(elements, 1)
    for alpha in 1:Ne
        inode = elements[alpha, :]
        x1 = nodes[inode[1], :]
        x2 = nodes[inode[2], :]
        x3 = nodes[inode[3], :]
        x4 = nodes[inode[4], :]
        V = tet_volume(x1, x2, x3, x4)

        for i in 1:4
            gi = 3*inode[i] - 2   # x 分量
            b[gi] += H0 * V / 4
        end
    end
end