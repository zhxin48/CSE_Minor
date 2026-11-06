using Test
using LinearAlgebra

# `@__DIR__` 是当前文件所在目录；`joinpath` 让路径同时适用于 Linux/Windows。
src = joinpath(@__DIR__, "..", "src")
include(joinpath(src, "volume.jl"))
include(joinpath(src, "mass.jl"))
include(joinpath(src, "load.jl"))
include(joinpath(src, "stiffness.jl"))
include(joinpath(src, "surface_stiffness.jl"))

const X1 = [0.0, 0.0, 0.0]
const X2 = [1.0, 0.0, 0.0]
const X3 = [0.0, 1.0, 0.0]
const X4 = [0.0, 0.0, 1.0]
const REF_TET = (X1, X2, X3, X4)
const FAR_TET = tuple((x .+ [5.0, 0.0, 0.0] for x in REF_TET)...)

@testset "geometry and local assembly" begin
    @test tet_volume(REF_TET...) ≈ 1 / 6
    @test size(basis_grads(REF_TET...)) == (3, 4)

    facets = @inferred tet_facets(REF_TET...)
    @test length(facets) == 4
    @test all(isapprox(norm(f[4]), 1.0; atol=1e-14) for f in facets)

    nodes = reduce(vcat, transpose.(REF_TET))
    elements = reshape([1, 2, 3, 4], 1, 4)
    A1 = zeros(12, 12)
    b = zeros(12)
    @test assemble_mass!(A1, elements, nodes, 2.0) === nothing
    @test assemble_load!(b, elements, nodes, 3.0) === nothing
    @test issymmetric(A1)
    @test sum(b) ≈ 3.0 / 6

    faces = @inferred boundary_faces(elements, nodes)
    @test length(faces) == 4
    @test all(dot(face[4], nodes[face[1], :] - nodes[setdiff(1:4, collect(face[1:3]))[1], :]) > 0
              for face in faces)
end

@testset "type-stable quadrature kernels" begin
    volume = @inferred integrate_tet(_ -> 1.0, REF_TET...;
                                     rtol=1e-3, maxevals=10_000)
    @test volume ≈ 1 / 6 atol=1e-12

    numerical = @inferred tet_tet_integral_num(
        REF_TET, FAR_TET; rtol=1e-2, atol=1e-8, maxevals=10_000)
    reduced = @inferred tet_tet_integral_facets(
        REF_TET, FAR_TET; rtol=1e-2, atol=1e-8, maxevals=10_000)
    @test numerical ≈ 0.005555531094601571 rtol=5e-4
    @test reduced ≈ 0.005555532350441779 rtol=5e-4
    @test numerical ≈ reduced rtol=5e-4

    self_matrix = @inferred self_facet_matrix(
        (X1, X2, X3); rtol=1e-2, atol=1e-8, maxevals=10_000)
    @test size(self_matrix) == (3, 3)
    @test all(isfinite, self_matrix)
    @test self_matrix ≈ self_matrix' rtol=5e-2
end

@testset "type-stable assembly path" begin
    nodes = reduce(vcat, transpose.(REF_TET))
    elements = reshape([1, 2, 3, 4], 1, 4)
    A2 = zeros(12, 12)

    volume_result = @inferred assemble_stiffness!(
        A2, elements, nodes;
        method=:facet, rtol=1e-2, atol=1e-8, maxevals=100_000)
    @test volume_result === A2

    surface_result = @inferred assemble_surface_stiffness!(
        A2, elements, nodes;
        rtol=1e-2, atol=1e-8, maxevals=100_000)
    @test surface_result === A2
    @test all(isfinite, A2)
end
