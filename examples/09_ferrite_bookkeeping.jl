# Ferrite prototype: replace hand-written FEM bookkeeping without changing the
# material or load formulas.
#
# Run from the repository root:
#   julia --project=. examples/09_ferrite_bookkeeping.jl

using LinearAlgebra

src = joinpath(@__DIR__, "..", "src")
include(joinpath(src, "volume.jl"))
include(joinpath(src, "mass.jl"))
include(joinpath(src, "load.jl"))
include(joinpath(src, "mesh.jl"))
include(joinpath(src, "ferrite_assembly.jl"))

nodes, elements = load_mesh()
chi_mag = 100.0
H_ext = (1.0, 0.0, 0.0)

# Ferrite now owns the grid topology, vector-valued P1 field, global degrees of
# freedom, cell iteration, quadrature data, and local-to-global assembly.
problem = ferrite_problem(nodes, elements)
A1_ferrite, b_ferrite = assemble_ferrite_material_load(
    problem; chi_mag=chi_mag, H_ext=H_ext)
boundary = ferrite_boundary_diagnostics(problem)

println("Ferrite grid:")
println("  cells           = ", getncells(problem.grid))
println("  dofs            = ", ndofs(problem.dofhandler))
println("  exterior facets = ", boundary.facet_count)
println("  exterior area   = ", boundary.area)
println("  integral(n dS)  = ", boundary.normal_integral)

# During migration the old implementation remains an independent reference.
ndof = 3 * size(nodes, 1)
A1_reference = zeros(ndof, ndof)
b_reference = zeros(ndof)
assemble_mass!(A1_reference, elements, nodes, chi_mag)
assemble_load!(b_reference, elements, nodes, H_ext[1])

println("\nAgreement with the existing implementation:")
println("  max row-sum norm of A1 difference = ",
        norm(Matrix(A1_ferrite) - A1_reference, Inf))
println("  max absolute b difference          = ",
        norm(b_ferrite - b_reference, Inf))


