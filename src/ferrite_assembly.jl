# Ferrite bookkeeping for mesh topology, local material terms, and the load vector.
#
# Ferrite 基础层：负责网格拓扑、向量自由度、局部材料项和载荷向量。
using Ferrite
using LinearAlgebra
using SparseArrays

# Build the Ferrite grid, topology, vector field, and public DOF permutation.
"""
    ferrite_problem(nodes, elements)

把仓库原有的 `nodes/elements` 数组转换为 Ferrite `Grid`，建立三分量 P1 磁化场，
并通过 topology 识别外边界。返回值同时保存 Ferrite 内部编号与仓库公开编号之间的映射。
"""
function ferrite_problem(nodes::AbstractMatrix{<:Real},
                         elements::AbstractMatrix{<:Integer})
    size(nodes, 2) == 3 || throw(ArgumentError("nodes must have three columns"))
    size(elements, 2) == 4 || throw(ArgumentError("elements must have four columns"))

    ferrite_nodes = [Node((Float64(nodes[i, 1]),
                           Float64(nodes[i, 2]),
                           Float64(nodes[i, 3]))) for i in axes(nodes, 1)]
    ferrite_cells = [Tetrahedron(Tuple(Int.(elements[e, :])))
                     for e in axes(elements, 1)]
    grid = Grid(ferrite_cells, ferrite_nodes)

    topology = ExclusiveTopology(grid)
    addboundaryfacetset!(grid, topology, "boundary", _ -> true)

    interpolation = Lagrange{RefTetrahedron, 1}()^3
    dofhandler = DofHandler(grid)
    add!(dofhandler, :M, interpolation)
    close!(dofhandler)

    # Ferrite 自主管理全局编号，不保证与输入节点顺序相同；因此显式保存到仓库公开顺序
    # `(Mx, My, Mz)` 的 permutation。
    nodal_dofs = zeros(Int, 3 * size(nodes, 1))
    for (cell_id, cell) in enumerate(ferrite_cells)
        cell_dofs = celldofs(dofhandler, cell_id)
        for local_node in 1:4, component in 1:3
            nodal_index = 3 * (cell.nodes[local_node] - 1) + component
            local_index = 3 * (local_node - 1) + component
            dof = cell_dofs[local_index]
            previous = nodal_dofs[nodal_index]
            (previous == 0 || previous == dof) || error("inconsistent Ferrite dof map")
            nodal_dofs[nodal_index] = dof
        end
    end
    all(!=(0), nodal_dofs) || error("some mesh nodes are not connected to a tetrahedron")

    public_dof = similar(nodal_dofs)
    for public_index in eachindex(nodal_dofs)
        public_dof[nodal_dofs[public_index]] = public_index
    end

    return (; grid, topology, dofhandler, interpolation, nodal_dofs, public_dof)
end

# Assemble the material mass matrix and external-field load with Ferrite.
"""
    assemble_ferrite_material_load(problem; chi_mag, H_ext)

使用 Ferrite 的 `CellValues`、`CellIterator`、`DofHandler` 和 `Assembler`
组装材料质量矩阵与外场载荷。返回结果转换回现有 `(Mx, My, Mz)` 节点交错顺序。
"""
function assemble_ferrite_material_load(problem;
                                        chi_mag::Real,
                                        H_ext=(1.0, 0.0, 0.0))
    chi_mag > 0 || throw(ArgumentError("chi_mag must be positive"))
    length(H_ext) == 3 || throw(ArgumentError("H_ext must have three components"))

    dh = problem.dofhandler
    quadrature = QuadratureRule{RefTetrahedron}(2)
    cellvalues = CellValues(quadrature, problem.interpolation)
    nlocal = getnbasefunctions(cellvalues)
    element_matrix = zeros(Float64, nlocal, nlocal)
    element_vector = zeros(Float64, nlocal)

    matrix = allocate_matrix(dh)
    vector = zeros(Float64, ndofs(dh))
    assembler = start_assemble(matrix, vector)
    applied_field = Vec{3, Float64}(Tuple(Float64.(H_ext)))

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(element_matrix, 0.0)
        fill!(element_vector, 0.0)

        for qp in 1:getnquadpoints(cellvalues)
            dV = getdetJdV(cellvalues, qp)
            for i in 1:nlocal
                Ni = shape_value(cellvalues, qp, i)
                element_vector[i] += dot(Ni, applied_field) * dV
                for j in 1:nlocal
                    Nj = shape_value(cellvalues, qp, j)
                    element_matrix[i, j] += dot(Ni, Nj) * dV / chi_mag
                end
            end
        end

        assemble!(assembler, celldofs(cell), element_matrix, element_vector)
    end

    # 转换回现有节点顺序，保证输出格式兼容旧代码和已有数据文件。
    p = problem.nodal_dofs
    return matrix[p, p], vector[p]
end

# Integrate boundary area and outward normals for topology diagnostics.
"""
    ferrite_boundary_diagnostics(problem)

在 Ferrite topology 得到的外边界 facet set 上积分面积与外法向。对闭合网格，
根据散度定理，`normal_integral` 应接近零。
"""
function ferrite_boundary_diagnostics(problem)
    facets = getfacetset(problem.grid, "boundary")
    quadrature = FacetQuadratureRule{RefTetrahedron}(2)
    facetvalues = FacetValues(quadrature, problem.interpolation)
    area = 0.0
    normal_integral = zeros(3)

    for facet in FacetIterator(problem.dofhandler, facets)
        reinit!(facetvalues, facet)
        for qp in 1:getnquadpoints(facetvalues)
            dS = getdetJdV(facetvalues, qp)
            area += dS
            normal_integral .+= getnormal(facetvalues, qp) .* dS
        end
    end

    return (; facet_count=length(facets), area, normal_integral)
end


