# Gmsh .msh reader: converts a four-node tetrahedral mesh into (nodes, elements).
using Gmsh

# Read a Gmsh .msh file and return (nodes, elements).
"""
    read_mesh(filename) -> (nodes, elements)

从 Gmsh `.msh` 文件读取网格 —— 既可用于外部 mesh，也可用于 `generate_box_mesh` 生成的标准 mesh。

- `nodes`：`Nn × 3` `Float64` 矩阵，节点坐标。
- `elements`：`Ne × 4` `Int` 矩阵，每行是一个四面体的 4 个局部节点索引。

Gmsh 的 node tag 不保证从 1 连续编号，因此读取时显式建立
`外部 tag → Julia 连续索引` 映射。当前只支持四节点线性四面体（Gmsh type 4）。
"""
function read_mesh(filename)
    Gmsh.initialize()
    try
        gmsh.open(filename)

        # 读取节点。Dict 推导式 `tag => i for ...` 建立外部标签到数组行号的映射。
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        node_index = Dict(tag => i for (i, tag) in enumerate(node_tags))
        Nn = length(node_tags)
        nodes = zeros(Float64, Nn, 3)
        for i in 1:Nn, component in 1:3
            nodes[i, component] = node_coords[3 * (i - 1) + component]
        end

        # `findfirst` 返回满足条件的第一个位置；Gmsh type 4 是四节点线性四面体。
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(3)
        tet_slot = findfirst(==(4), elem_types)
        tet_slot === nothing && error("mesh 中没有四节点线性四面体（Gmsh element type 4）")

        Ne = length(elem_tags[tet_slot])
        elements = zeros(Int, Ne, 4)
        connectivity = elem_node_tags[tet_slot]
        for element in 1:Ne, local_node in 1:4
            tag = connectivity[4 * (element - 1) + local_node]
            elements[element, local_node] = node_index[tag]
        end

        return nodes, elements
    finally
        # 即使读取或格式检查失败，也要释放 Gmsh 的全局状态。
        Gmsh.finalize()
    end
end
