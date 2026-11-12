# 主入口：加载 mesh → 组装 A1/A2/A、b → 求解 A·u = b → 物理解验 → 保存 + 可视化。
#
# 用法：
#   julia --project=. src/main.jl                 # 标准 mesh（data/my_box.msh，不存在则生成）
#   julia --project=. src/main.jl path/to.msh     # 外部 Gmsh mesh
include("mesh.jl")       # load_mesh（含 generate / read）
include("volume.jl")
include("mass.jl")
include("load.jl")
include("stiffness.jl")
include("surface_stiffness.jl")
include("solve.jl")
include("visualize.jl")

function main()
    # 1. mesh：无参数用标准 mesh，有参数加载外部 .msh
    filename = isempty(ARGS) ? nothing : ARGS[1]
    nodes, elements = load_mesh(filename=filename)
    Nn = size(nodes, 1)
    Ne = size(elements, 1)
    ndof = 3 * Nn
    println("使用 mesh: ", filename === nothing ? "标准 mesh" : filename)
    println("Nn=$Nn  Ne=$Ne  ndof=$ndof")

    chi_mag = 100.0
    H0 = 1.0

    # 2. 组装
    A1 = zeros(ndof, ndof); assemble_mass!(A1, elements, nodes, chi_mag)
    A2 = zeros(ndof, ndof); assemble_stiffness!(A2, elements, nodes; precision=:fast)  # VV 体积散度项
    # 表面磁荷项 SS − SV − VS（自面积分用 Duffy 变换）。立方体 trace=1、单四面体 Nxx=1/3 已验证。
    # 注意：标准 1×1×0.1 板网格只有 2 层 z 单元，退磁因子偏小（mean Mx≈12 vs 解析 9.3），属网格分辨率。
    assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-3)                         # + SS − SV − VS
    A = A1 + A2
    b = zeros(ndof);        assemble_load!(b, elements, nodes, H0)
    println("sum(A2) = ", sum(A2), "（含表面磁荷项后不再 ≈ 0）")

    # 3. 求解
    u = solve_system(A, b)

    # 4. 物理解验
    verify_solution(A, b, u, nodes, elements; H0=H0, chi_mag=chi_mag)

    # 5. 保存 + 可视化
    datadir = joinpath(@__DIR__, "..", "data")
    plotdir = joinpath(@__DIR__, "..", "plots")
    save_solution(u, joinpath(datadir, "u.txt"))
    plot_mesh(nodes, elements; out=joinpath(plotdir, "mesh.png"))
    plot_field(nodes, elements, u; out=joinpath(plotdir, "field.png"))
    plot_field_components(nodes, elements, u; outdir=plotdir)
end

main()
