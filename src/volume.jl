using LinearAlgebra

"""
    tet_volume(x1, x2, x3, x4)

计算四面体体积。
x1, x2, x3, x4 是长度为 3 的向量，表示四个顶点坐标。
"""
function tet_volume(x1, x2, x3, x4)
    v1 = x2 - x1
    v2 = x3 - x1
    v3 = x4 - x1
    M = hcat(v1, v2, v3)   # 3×3 矩阵，每列是一个向量
    return abs(det(M)) / 6
end