# 统一 mesh 入口：把「标准 mesh 生成」和「外部 mesh 加载」收敛到一个函数 load_mesh。
include("generate_mesh.jl")
include("read_mesh.jl")

const STANDARD_MESH_PATH = joinpath(@__DIR__, "..", "data", "my_box.msh")

"""
    load_mesh(; filename=nothing, force=false, kwargs...) -> (nodes, elements)

统一网格入口：

- 给定 `filename`（外部 Gmsh `.msh` 路径）→ `read_mesh(filename)`。
- 未给定 → 使用标准 mesh `data/my_box.msh`：文件不存在或 `force=true` 时先用
  `generate_box_mesh` 生成（`kwargs` 传入生成器，如 `lx/ly/lz/mesh_size_min/mesh_size_max`），再读取。

用法：
    nodes, elements = load_mesh()                              # 标准 mesh
    nodes, elements = load_mesh(filename="my/external.msh")    # 外部 mesh
    nodes, elements = load_mesh(force=true, lx=2.0, lz=0.2)    # 重新生成标准 mesh（改尺寸）
"""
function load_mesh(; filename=nothing, force=false, kwargs...)
    # `;` 后是关键字参数；`kwargs...` 收集其余命名参数并转交网格生成器。
    # `condition ? a : b` 是三元表达式：条件成立取 a，否则取 b。
    path = filename === nothing ? STANDARD_MESH_PATH : filename
    if filename === nothing && (force || !isfile(path))
        generate_box_mesh(path; kwargs...)
    end
    return read_mesh(path)
end
