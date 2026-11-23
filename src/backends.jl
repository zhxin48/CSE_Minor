# Dual bookkeeping/assembly backends for Part B.
#
# Part B 双后端：两条路径共享 HCubature、Euler 和 Duffy 数学内核，只改变 cells、
# facets、coordinates 与 DOFs 的来源和缓存方式。
using LinearAlgebra
using Ferrite

abstract type AssemblyBackend end

struct LegacyBackend{TN, TE, TC, TF} <: AssemblyBackend
    nodes::TN
    elements::TE
    cells::TC
    boundary::TF
end

struct FerriteBackend{TN, TE, TP, TC, TF} <: AssemblyBackend
    nodes::TN
    elements::TE
    problem::TP
    cells::TC
    boundary::TF
end

# Return the repository's three public DOFs for one node.
# 返回一个节点在仓库公开交错顺序中的三个自由度。
_public_node_dofs(node::Integer) = ntuple(component -> 3 * (node - 1) + component, 3)

# Build cell records from the hand-written arrays.
# 从手写 `nodes/elements` 数组构造体单元记录。
function _legacy_cell_records(nodes, elements)
    return [begin
        node_ids = Tuple(Int.(elements[cell_id, :]))
        coordinates = ntuple(i -> Vector{Float64}(nodes[node_ids[i], :]), 4)
        dofs = [dof for node in node_ids for dof in _public_node_dofs(node)]
        (; id=cell_id, nodes=node_ids, coordinates, dofs)
    end for cell_id in axes(elements, 1)]
end

# Build deterministic boundary records with the legacy face counter.
# 使用旧面计数器构造确定性排序的边界记录。
function _legacy_boundary_records(nodes, elements)
    records = [begin
        i, j, k, normal = face
        node_ids = (i, j, k)
        coordinates = ntuple(q -> Vector{Float64}(nodes[node_ids[q], :]), 3)
        dofs = [dof for node in node_ids for dof in _public_node_dofs(node)]
        (; nodes=node_ids, coordinates, normal=Vector{Float64}(normal), dofs)
    end for face in boundary_faces(elements, nodes)]
    sort!(records; by=record -> Tuple(sort(collect(record.nodes))))
    return records
end

# Construct the complete legacy bookkeeping backend.
# 构造完整手写 bookkeeping 后端。
function LegacyBackend(nodes::AbstractMatrix{<:Real},
                       elements::AbstractMatrix{<:Integer})
    cells = _legacy_cell_records(nodes, elements)
    boundary = _legacy_boundary_records(nodes, elements)
    return LegacyBackend(nodes, elements, cells, boundary)
end

# Build concrete cell records through Ferrite CellIterator.
# 通过 Ferrite `CellIterator` 构造具体类型的体单元记录。
function _ferrite_cell_records(problem)
    records = [begin
        node_ids = Tuple(Int.(getnodes(cell)))
        coordinates = Tuple(Vector{Float64}(coordinate) for coordinate in getcoordinates(cell))
        dofs = problem.public_dof[collect(celldofs(cell))]
        (; id=cellid(cell), nodes=node_ids, coordinates, dofs)
    end for cell in CellIterator(problem.dofhandler)]
    sort!(records; by=record -> record.id)
    return records
end

# Build deterministic exterior-facet records through Ferrite topology.
# 通过 Ferrite topology 构造确定性排序的外边界 facet 记录。
function _ferrite_boundary_records(problem)
    grid = problem.grid
    facetvalues = FacetValues(FacetQuadratureRule{RefTetrahedron}(1),
                              problem.interpolation)
    cache = FacetCache(problem.dofhandler)
    records = [begin
        cell_id, local_facet = facet_index
        reinit!(cache, facet_index)
        reinit!(facetvalues, cache)

        cell = getcells(grid, cell_id)
        node_ids = Tuple(Int.(Ferrite.facets(cell)[local_facet]))
        coordinates = ntuple(i -> Vector{Float64}(grid.nodes[node_ids[i]].x), 3)
        normal = Vector{Float64}(getnormal(facetvalues, 1))

        cell_nodes = collect(getnodes(cache))
        cell_dofs = problem.public_dof[collect(celldofs(cache))]
        dofs = Int[]
        for node in node_ids
            local_node = findfirst(==(node), cell_nodes)
            local_node === nothing && error("facet node is not part of its owner cell")
            append!(dofs, cell_dofs[(3local_node - 2):(3local_node)])
        end
        (; nodes=node_ids, coordinates, normal, dofs,
           owner=cell_id, local_facet)
    end for facet_index in getfacetset(grid, "boundary")]
    sort!(records; by=record -> Tuple(sort(collect(record.nodes))))
    return records
end

# Construct the complete Ferrite bookkeeping backend.
# 构造完整 Ferrite bookkeeping 后端。
function FerriteBackend(nodes::AbstractMatrix{<:Real},
                        elements::AbstractMatrix{<:Integer})
    problem = ferrite_problem(nodes, elements)
    cells = _ferrite_cell_records(problem)
    boundary = _ferrite_boundary_records(problem)
    return FerriteBackend(nodes, elements, problem, cells, boundary)
end

# Return the user-facing backend name and problem size.
# 返回用于输出的后端名称和问题自由度数。
backend_name(::LegacyBackend) = "legacy"
backend_name(::FerriteBackend) = "ferrite"
backend_ndofs(backend::AssemblyBackend) = 3 * size(backend.nodes, 1)

# Assemble material and load terms with the hand-written path.
# 使用手写路径组装材料质量矩阵和载荷向量。
function assemble_backend_material_load(backend::LegacyBackend;
                                        chi_mag::Real,
                                        H_ext=(1.0, 0.0, 0.0))
    H_ext[2] == 0 && H_ext[3] == 0 ||
        throw(ArgumentError("legacy load assembly currently supports an x-directed field"))
    A1 = zeros(backend_ndofs(backend), backend_ndofs(backend))
    b = zeros(backend_ndofs(backend))
    assemble_mass!(A1, backend.elements, backend.nodes, chi_mag)
    assemble_load!(b, backend.elements, backend.nodes, H_ext[1])
    return A1, b
end

# Assemble material and load terms with Ferrite CellValues and Assembler.
# 使用 Ferrite `CellValues` 与 `Assembler` 组装材料项和载荷项。
function assemble_backend_material_load(backend::FerriteBackend;
                                        chi_mag::Real,
                                        H_ext=(1.0, 0.0, 0.0))
    A1, b = assemble_ferrite_material_load(
        backend.problem; chi_mag=chi_mag, H_ext=H_ext)
    return Matrix(A1), b
end

# Select the numerical or Euler-reduced kernel for one cell pair.
# 为一个 cell pair 选择直接数值积分或 Euler 降维积分。
function _use_facet_method(method, alpha, beta, cell_a, cell_b)
    method == :facet && return true
    method == :num && return false
    method == :hybrid && return alpha == beta
    method == :hybrid2 && return alpha == beta || !isempty(intersect(cell_a.nodes, cell_b.nodes))
    error("unknown method: $method")
end

# Evaluate one volume-volume interaction block with shared mathematical kernels.
# 使用共享数学内核计算一个体积--体积局部块。
function _volume_pair_block(cell_a, cell_b;
                            method, rtol, atol, maxevals)
    use_facet = _use_facet_method(method, cell_a.id, cell_b.id, cell_a, cell_b)
    integral = if use_facet
        tet_tet_integral_facets(cell_a.coordinates, cell_b.coordinates;
                                rtol=rtol, atol=atol, maxevals=maxevals)
    else
        tet_tet_integral_num(cell_a.coordinates, cell_b.coordinates;
                             rtol=rtol, atol=atol, maxevals=maxevals)
    end
    gradients_a = basis_grads(cell_a.coordinates...)
    gradients_b = basis_grads(cell_b.coordinates...)
    block = zeros(12, 12)
    for i in 1:4, j in 1:4, a in 1:3, b in 1:3
        block[3(i - 1) + a, 3(j - 1) + b] =
            gradients_a[a, i] * gradients_b[b, j] * integral / (4π)
    end
    return block
end

# Assemble the complete VV block from backend-provided cell records.
# 从后端提供的 cell records 组装完整 `VV` 块。
function assemble_backend_vv(backend::AssemblyBackend;
                             precision=:fast, method=nothing,
                             rtol=nothing, atol=nothing,
                             maxevals::Integer=10^7,
                             progress::Bool=false)
    preset_method, preset_rtol, preset_atol = _stiffness_preset(precision)
    method = something(method, preset_method)
    rtol = something(rtol, preset_rtol)
    atol = something(atol, preset_atol)
    method == :num && error("method=:num cannot assemble self interactions")

    A = zeros(backend_ndofs(backend), backend_ndofs(backend))
    ncells = length(backend.cells)
    total = ncells * (ncells + 1) ÷ 2
    iteration = 0
    for alpha in 1:ncells, beta in alpha:ncells
        iteration += 1
        progress && _progress(iteration, total, "$(backend_name(backend)) VV")
        cell_a = backend.cells[alpha]
        cell_b = backend.cells[beta]
        block = _volume_pair_block(cell_a, cell_b;
                                   method=method, rtol=rtol, atol=atol,
                                   maxevals=maxevals)
        A[cell_a.dofs, cell_b.dofs] .+= block
        alpha == beta || (A[cell_b.dofs, cell_a.dofs] .+= block')
    end
    return A
end

# Assemble SS and the mixed SV/VS blocks from backend facet records.
# 从后端提供的 facet records 组装 `SS` 与混合 `SV/VS` 块。
function assemble_backend_surface_blocks(backend::AssemblyBackend;
                                         rtol::Real=1e-4,
                                         atol::Real=1e-8,
                                         maxevals::Integer=10^7,
                                         progress::Bool=false)
    ndof = backend_ndofs(backend)
    SS = zeros(ndof, ndof)
    MIXED = zeros(ndof, ndof) # -(SV + VS)
    faces = backend.boundary
    cells = backend.cells
    inv4pi = 1 / (4π)

    total_ss = length(faces) * (length(faces) + 1) ÷ 2
    iteration = 0
    for a in eachindex(faces), b in a:length(faces)
        iteration += 1
        progress && _progress(iteration, total_ss, "$(backend_name(backend)) SS")
        face_a = faces[a]
        face_b = faces[b]
        S = a == b ?
            self_facet_matrix(face_a.coordinates; rtol=rtol, atol=atol, maxevals=maxevals) :
            facet_facet_matrix(face_a.coordinates, face_b.coordinates;
                               rtol=rtol, atol=atol, maxevals=maxevals)
        a == b && (S = (S + S') / 2)
        for m in 1:3, n in 1:3, component_a in 1:3, component_b in 1:3
            row = face_a.dofs[3(m - 1) + component_a]
            col = face_b.dofs[3(n - 1) + component_b]
            value = inv4pi * face_a.normal[component_a] *
                    face_b.normal[component_b] * S[m, n]
            SS[row, col] += value
            a == b || (SS[col, row] += value)
        end
    end

    total_mixed = length(faces) * length(cells)
    iteration = 0
    cell_facets = [tet_facets(cell.coordinates...) for cell in cells]
    cell_gradients = [basis_grads(cell.coordinates...) for cell in cells]
    for face in faces, alpha in eachindex(cells)
        iteration += 1
        progress && _progress(iteration, total_mixed, "$(backend_name(backend)) SV/VS")
        V = zeros(3)
        for (p1, p2, p3, normal) in cell_facets[alpha]
            V .+= facet_facet_sv_vector(face.coordinates, (p1, p2, p3), normal;
                                        rtol=rtol, atol=atol, maxevals=maxevals)
        end
        V .*= 0.5

        cell = cells[alpha]
        gradients = cell_gradients[alpha]
        for m in 1:3, local_node in 1:4, component_face in 1:3, component_cell in 1:3
            value = inv4pi * gradients[component_cell, local_node] *
                    face.normal[component_face] * V[m]
            cell_dof = cell.dofs[3(local_node - 1) + component_cell]
            face_dof = face.dofs[3(m - 1) + component_face]
            MIXED[cell_dof, face_dof] -= value
            MIXED[face_dof, cell_dof] -= value
        end
    end
    return (; SS, MIXED, total=SS + MIXED)
end

# Assemble all four contributions of the demagnetising operator.
# 组装完整四项退磁算子。
function assemble_backend_demag(backend::AssemblyBackend;
                                precision=:fast, method=nothing,
                                vv_rtol=nothing, vv_atol=nothing,
                                surface_rtol::Real=1e-4,
                                surface_atol::Real=1e-8,
                                maxevals::Integer=10^7,
                                progress::Bool=false)
    VV = assemble_backend_vv(backend; precision, method,
                             rtol=vv_rtol, atol=vv_atol,
                             maxevals, progress)
    surface = assemble_backend_surface_blocks(backend;
                                              rtol=surface_rtol,
                                              atol=surface_atol,
                                              maxevals, progress)
    return (; VV, SS=surface.SS, MIXED=surface.MIXED,
             total=VV + surface.total)
end

# Assemble and solve the complete linear system for one backend.
# 通过指定后端组装并求解完整线性系统。
function solve_backend(backend::AssemblyBackend;
                       chi_mag::Real=100.0,
                       H_ext=(1.0, 0.0, 0.0), kwargs...)
    A1, b = assemble_backend_material_load(backend; chi_mag, H_ext)
    demag = assemble_backend_demag(backend; kwargs...)
    A = A1 + demag.total
    u = A \ b
    return (; A1, b, demag, A, u)
end

