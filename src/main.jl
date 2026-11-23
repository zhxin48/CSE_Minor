# Main entry point: choose a bookkeeping backend, assemble, solve, verify, and plot.
#
# 主入口：选择 bookkeeping 后端，加载网格，组装完整系统，求解、验证并输出图片。
# 默认使用 Ferrite；`--backend=legacy` 完整保留原来的手写装配路径。
#
# 用法：
#   julia --project=. src/main.jl
#   julia --project=. src/main.jl --backend=ferrite
#   julia --project=. src/main.jl --backend=legacy
#   julia --project=. src/main.jl --backend=ferrite path/to/mesh.msh
include("mesh.jl")
include("volume.jl")
include("mass.jl")
include("load.jl")
include("stiffness.jl")
include("surface_stiffness.jl")
include("ferrite_assembly.jl")
include("backends.jl")
include("solve.jl")
include("visualize.jl")

# Parse the backend option and the optional external mesh path.
"""
    parse_main_args(args) -> (backend, filename)

解析命令行参数。backend 可选 ferrite 或 legacy；没有 mesh 路径时使用标准网格。
同时支持 --backend=ferrite 和 --backend ferrite 两种写法。
"""
function parse_main_args(args)
    backend = :ferrite
    filename = nothing
    index = 1
    while index <= length(args)
        argument = args[index]
        if startswith(argument, "--backend=")
            backend = Symbol(lowercase(split(argument, "="; limit=2)[2]))
        elseif argument == "--backend"
            index == length(args) && error("--backend 后需要 ferrite 或 legacy")
            index += 1
            backend = Symbol(lowercase(args[index]))
        elseif startswith(argument, "--")
            error("未知选项：$argument")
        elseif filename === nothing
            filename = argument
        else
            error("只能提供一个外部 mesh 路径")
        end
        index += 1
    end
    backend in (:ferrite, :legacy) || error("backend 必须是 ferrite 或 legacy")
    return backend, filename
end

# Assemble the complete system with the original hand-written bookkeeping.
"""
    assemble_legacy_system(nodes, elements, chi_mag, H0)

使用迁移前的完整手写路径：手动 DOF 编号、boundary_faces、体积项和表面项装配。
该函数作为可随时回退的稳定基线，不经过 Ferrite 数据结构。
"""
function assemble_legacy_system(nodes, elements, chi_mag, H0)
    ndof = 3 * size(nodes, 1)
    A1 = zeros(ndof, ndof)
    assemble_mass!(A1, elements, nodes, chi_mag)

    A2 = zeros(ndof, ndof)
    assemble_stiffness!(A2, elements, nodes; precision=:fast)
    assemble_surface_stiffness!(A2, elements, nodes; rtol=1e-3)

    b = zeros(ndof)
    assemble_load!(b, elements, nodes, H0)
    return A1, A2, b
end

# Assemble the complete system with Ferrite topology and degree-of-freedom data.
"""
    assemble_ferrite_system(nodes, elements, chi_mag, H0)

使用 Ferrite 管理 Grid、topology、DOF、cells、facets 和 scatter；非局部积分仍调用
项目已经验证的 HCubature、Euler 和 Duffy 数学内核。
"""
function assemble_ferrite_system(nodes, elements, chi_mag, H0)
    backend = FerriteBackend(nodes, elements)
    A1, b = assemble_backend_material_load(
        backend; chi_mag=chi_mag, H_ext=(H0, 0.0, 0.0))
    demag = assemble_backend_demag(
        backend; precision=:fast, surface_rtol=1e-3,
        maxevals=10^7, progress=true)
    return A1, demag.total, b
end

# Run the selected complete solver pipeline.
"""
    main(args=ARGS)

运行完整求解流程。默认后端为 Ferrite；选择 legacy 时仍执行迁移前的手写实现。
"""
function main(args=ARGS)
    backend, filename = parse_main_args(args)
    nodes, elements = load_mesh(filename=filename)
    Nn = size(nodes, 1)
    Ne = size(elements, 1)
    ndof = 3 * Nn
    println("使用 mesh: ", filename === nothing ? "标准 mesh" : filename)
    println("backend=$(backend)  Nn=$Nn  Ne=$Ne  ndof=$ndof")

    chi_mag = 100.0
    H0 = 1.0
    A1, A2, b = backend == :legacy ?
        assemble_legacy_system(nodes, elements, chi_mag, H0) :
        assemble_ferrite_system(nodes, elements, chi_mag, H0)
    A = A1 + A2
    println("sum(A2) = ", sum(A2), "（完整退磁算子）")

    u = solve_system(A, b)
    verify_solution(A, b, u, nodes, elements; H0=H0, chi_mag=chi_mag)

    datadir = joinpath(@__DIR__, "..", "data")
    plotdir = joinpath(@__DIR__, "..", "plots")
    save_solution(u, joinpath(datadir, "u.txt"))
    plot_mesh(nodes, elements; out=joinpath(plotdir, "mesh.png"))
    plot_field(nodes, elements, u; out=joinpath(plotdir, "field.png"))
    plot_field_components(nodes, elements, u; outdir=plotdir)
    return u
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
