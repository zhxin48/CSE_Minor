# Magnetized Metallic Plate

Run all commands from the project root. Only `Project.toml` is shared by Git;
each user resolves a local `Manifest.toml`, which is intentionally ignored.

## Install the Julia environment

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Resolves versions allowed by the compatibility bounds in `Project.toml`, creates
a local `Manifest.toml`, and precompiles the packages.

If an old local environment cannot be resolved, delete only the local manifest
and run the command again:

```bash
rm -f Manifest.toml
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

## Use the interactive Julia terminal

Start Julia from the WSL shell with this project active:

```bash
cd /home/xin/CSE_Minor
julia --project=.
```

Then enter Julia commands at the `julia>` prompt:

```julia
using Pkg
Pkg.instantiate()
Pkg.precompile()
include("src/main.jl")
main(["--backend=ferrite"])
```

Use `main(["--backend=legacy"])` for the complete hand-written fallback, or
`main(["--backend=ferrite", "path/to/mesh.msh"])` for an external mesh. Run
`exit()` or press `Ctrl-D` to leave Julia. Shell commands such as
`julia --project=.` must be entered before the `julia>` prompt, not inside it.

## Run the tests

```bash
julia --project=. test/runtests.jl
```

Tests geometry, assembly, type stability, HCubature, Euler reduction, and Duffy integration.

## Generate the standard mesh

```bash
julia --project=. src/generate_mesh.jl
```

Generates the default Gmsh tetrahedral mesh at `data/my_box.msh`.

## Run the complete solver

```bash
julia --project=. src/main.jl
```

Uses the Ferrite bookkeeping backend by default, assembles and solves the full
system, then writes `data/u.txt` and images to `plots/`.

The Ferrite backend manages mesh topology, vector degrees of freedom, cells,
exterior facets, outward normals, and local-to-global mappings. The established
HCubature, Euler, and Duffy integration kernels remain project-specific.

## Select the bookkeeping backend

```bash
# Default: Ferrite Grid, topology, DofHandler, CellIterator, and FacetIterator
julia --project=. src/main.jl --backend=ferrite

# Complete pre-migration path: hand-written DOFs, boundary_faces, and assembly
julia --project=. src/main.jl --backend=legacy
```

`legacy` is a genuine fallback: it calls the original `assemble_mass!`,
`assemble_load!`, `assemble_stiffness!`, and `assemble_surface_stiffness!`
functions without constructing Ferrite data structures.

## Solve with an external mesh

```bash
julia --project=. src/main.jl --backend=ferrite path/to/mesh.msh
julia --project=. src/main.jl --backend=legacy path/to/mesh.msh
```

Runs the complete solver using a four-node tetrahedral Gmsh mesh.

## Inspect the Ferrite bookkeeping layer

```bash
julia --project=. examples/09_ferrite_bookkeeping.jl
```

Prints the Ferrite cell, DOF, exterior-facet, area, and closed-surface normal
diagnostics and compares the Ferrite material/load assembly with the original
implementation.

## Solve on a structured mesh

```bash
julia --project=. src/scripts/solve_structured.jl
```

Runs the same physical model on a regular mesh with explicitly controlled thickness layers.

## Generate and save the matrices

```bash
julia --project=. src/scripts/build_matrices.jl
```

Writes `data/A1.txt`, `data/A2.txt`, and `data/A.txt`.

## Check the volume stiffness assembly

```bash
julia --project=. src/scripts/assemble_global.jl
```

Assembles the `VV` contribution and reports symmetry, row sums, and matrix statistics.

## Verify the integration methods

```bash
julia --project=. src/scripts/check_stiffness.jl
```

Compares methods A and B and checks the tetrahedron self-interaction.

## Run the performance benchmark

```bash
julia --project=. src/scripts/benchmark.jl
```

Compares the accuracy and runtime of the numerical and dimension-reduced methods.

## Verify the Aharoni demagnetizing factor

```bash
julia --project=. src/scripts/converge_thickness.jl
```

Checks whether the discrete surface formulation converges to the Aharoni value of \(N_x\).

## Verify convergence of the full solution

```bash
julia --project=. src/scripts/verify_mean_mx.jl
```

Checks whether the mean \(M_x\) of the non-uniform self-consistent solution converges under mesh refinement.

## Run the detailed diagnostics

```bash
julia --project=. src/scripts/dlahaye_followup_checks.jl
```

Checks type inference, nested HCubature costs, boundary faces, and field-component plots.
