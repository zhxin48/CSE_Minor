using Gmsh

"""
    generate_box_mesh(filename; lx=1.0, ly=1.0, lz=0.1, mesh_size_min=0.2, mesh_size_max=0.5, algorithm3d=1)

生成一个长方体 `[0,lx]×[0,ly]×[0,lz]` 的四面体网格，写入 Gmsh `.msh` 文件。

这是「标准 mesh」的生成器。要加载外部 mesh，用 `read_mesh`（或 `load_mesh(; filename=...)`）。
"""
function generate_box_mesh(filename;
                           lx=1.0, ly=1.0, lz=0.1,
                           mesh_size_min=0.2, mesh_size_max=0.5,
                           algorithm3d=1)
    d = dirname(filename)
    isempty(d) || mkpath(d)

    Gmsh.initialize()

    # 网格尺寸与 3D 算法（1 = Delaunay）
    gmsh.option.setNumber("General.Verbosity", 3)
    gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)

    # 几何：长方体
    gmsh.model.occ.add_box(0, 0, 0, lx, ly, lz, 1)
    gmsh.model.occ.synchronize()

    # 生成四面体网格
    gmsh.model.mesh.generate(3)

    gmsh.write(filename)
    Gmsh.finalize()

    println("Mesh written to: ", filename)
end

# 直接运行：julia --project=. src/generate_mesh.jl [输出路径]
if abspath(PROGRAM_FILE) == @__FILE__
    out = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data", "my_box.msh")
    generate_box_mesh(out)
end
