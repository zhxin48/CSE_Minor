# Full pipeline (assemble → solve → physical checks) on a structured mesh.
# Usage: julia --project=. src/scripts/solve_structured.jl
# Result: mean Mx ≈ 12.2 (non-uniform solution, physically correct; 9.3 is only
# the uniform-magnetisation approximation, not the exact solution).
#
# 用结构化网格做全流程（组装→求解→物理检验）。
# 用法：julia --project=. src/scripts/solve_structured.jl
# 结果：mean Mx ≈ 12.2（非均匀解，物理正确；9.3 只是均匀磁化近似，不是精确解）。
include(joinpath(@__DIR__, "..", "structured_mesh.jl"))
include(joinpath(@__DIR__, "..", "volume.jl"))
include(joinpath(@__DIR__, "..", "mass.jl"))
include(joinpath(@__DIR__, "..", "load.jl"))
include(joinpath(@__DIR__, "..", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "surface_stiffness.jl"))
include(joinpath(@__DIR__, "..", "solve.jl"))
using LinearAlgebra

# Run the full pipeline on an nx×ny×nz structured mesh.
function main()
    nx, ny, nz = 6, 6, 2          # 结构化网格（nx×ny×nz 个单元）
    nodes, elements = structured_box_mesh(nx, ny, nz)
    Nn = size(nodes, 1)
    ndof = 3 * Nn
    println("结构化网格：nx=", nx, " ny=", ny, " nz=", nz, "  Nn=", Nn, "  Ne=", size(elements, 1))

    chi_mag = 100.0
    H0 = 1.0

    A1 = zeros(ndof, ndof); assemble_mass!(A1, elements, nodes, chi_mag)
    A2 = zeros(ndof, ndof); assemble_stiffness!(A2, elements, nodes; precision=:fast)   # VV 体积项
    assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-3)                          # SS − SV − VS 表面项
    A = A1 + A2
    b = zeros(ndof);        assemble_load!(b, elements, nodes, H0)

    u = solve_system(A, b)
    verify_solution(A, b, u, nodes, elements; H0=H0, chi_mag=chi_mag)
end

main()
