# 专题 2：结构化网格、Gmsh 网格和边界面。

include(joinpath(@__DIR__, "..", "src", "structured_mesh.jl"))
include(joinpath(@__DIR__, "..", "src", "mesh.jl"))
include(joinpath(@__DIR__, "..", "src", "stiffness.jl"))
include(joinpath(@__DIR__, "..", "src", "surface_stiffness.jl"))

println("=== A. 结构化网格 ===")
structured_nodes, structured_elements = structured_box_mesh(
    2, 2, 1; lx=1.0, ly=1.0, lz=0.1)
structured_boundary = boundary_faces(structured_elements, structured_nodes)
println("nodes = ", size(structured_nodes, 1))
println("tetrahedra = ", size(structured_elements, 1))
println("boundary triangles = ", length(structured_boundary))
println("first node = ", structured_nodes[1, :])
println("first tetrahedron = ", structured_elements[1, :])

println("\n=== B. Gmsh 标准网格 ===")
gmsh_nodes, gmsh_elements = load_mesh()
gmsh_boundary = boundary_faces(gmsh_elements, gmsh_nodes)
println("nodes = ", size(gmsh_nodes, 1))
println("tetrahedra = ", size(gmsh_elements, 1))
println("boundary triangles = ", length(gmsh_boundary))
println("local index range = ", extrema(gmsh_elements))

println("""

read_mesh 的职责：
1. 读取 Gmsh 的外部 node tag 和坐标；
2. 建立 tag => 1:Nn 的 Dict；
3. 把单元连接关系改写成 Julia 连续局部索引；
4. 只把 nodes/elements 交给后续数值模块。

因此后续装配代码不需要知道原始 tag 是否连续。
""")
