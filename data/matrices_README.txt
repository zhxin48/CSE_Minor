本目录矩阵由 src/scripts/build_matrices.jl 生成。
运行：julia --project=. src/scripts/build_matrices.jl

网格：my_box.msh（Nn=108, Ne=262，ndof=324）
物理参数：chi_mag=100.0, H_ext=(H0,0,0), H0=1.0
积分精度：precision=fast（method=hybrid, rtol=1e-3, atol=1e-6）

A1.txt  — 质量矩阵 A^(1) = (1/χ_mag)·∫ φ_j·φ_k dΩ（分块对角）
A2.txt  — 刚度矩阵 A^(2) = (1/4π)[VV − SV − VS + SS]（体积散度 + 表面磁荷，完整对称弱形式）
A.txt   — 合并 A = A1 + A2，线性系统 A·u = b 的系数矩阵

读回 Julia：using DelimitedFiles; A = readdlm("A.txt")
