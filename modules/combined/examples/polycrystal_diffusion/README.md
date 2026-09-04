# Polycrystal grain-boundary diffusion

This directory is a runnable proof of concept for diffusion through a Neper polycrystal. Grain
interiors are ordinary finite elements and grain boundaries are conforming lower-dimensional
elements. Both use the same nodal concentration `c`; there is no GB concentration, exchange
coefficient, mortar transfer, or geometrically resolved GB thickness.

## Model and units

The 3D weak form is

\[
\int_\Omega D_b\nabla c\cdot\nabla v\,dV+
\int_{\Gamma_{GB}}K_{GB}^{ex}\nabla_\Gamma c\cdot\nabla_\Gamma v\,dA=0,
\]

with the analogous area-plus-line form in 2D. The supplied demonstration values are

```text
D_bulk       = 1e-19 m^2/s
D_gb         = 1e-13 m^2/s = 1e6 D_bulk
gb_thickness = 0.5 nm
K_gb_excess  = gb_thickness * (D_gb - D_bulk)
             = 4.999995e-23 m^3/s
```

`K_gb_excess`, rather than `gb_thickness*D_gb`, avoids counting the bulk contribution twice.
MOOSE `${units ...}` expressions document and convert every dimensional input to SI. Neper still
works in micrometers for robust meshing; `TransformGenerator` scales its mesh by `1e-6` before the
solve.

`MatDiffusion` supplies the bulk and GB terms. On an embedded EDGE2 or TRI3, libMesh maps the
shape-function gradient into the tangent line or plane, so the second term is intrinsic GB
diffusion. `TimeDerivative` is restricted to the bulk because this model adds conductance but no
excess GB storage. All executioners use MOOSE automatic scaling.

## Finalized inputs

- `polycrystal_diffusion.i`: steady 2D flux verification.
- `polycrystal_transient.i`: 2D, 100-grain constant-source transient through 200 um by 50 um.
- `bicrystal_whipple.i`: fine 2D bicrystal checked against Whipple's constant-source solution.
- `polycrystal_3d.i`: quadratic 3D, 100-grain transient through 200 um by 50 um by 50 um. The complete
  square `x=0` face has `c=1`, the complete square `x=200 um` face has `c=0`, and the other four
  faces have natural zero flux.

The default 2D mesh contains TRI3 bulk and EDGE2 GB elements. Neper generates TET4/TRI3 in 3D,
then MOOSE promotes the complete mixed-dimensional mesh to TET10 bulk and TRI6 GB elements. This
substantially suppresses the early-time undershoot on the coarse tetrahedral topology.
`LowerDBlockFromSidesetGenerator` converts only the internal `grain_boundaries` sideset to a
lower-dimensional block and reuses the parent-element nodes.

## Run

From this directory:

```bash
# 2D transient, 1e9 s
conda run -n moose ../../../../test/moose_test-opt -i polycrystal_transient.i

# 2D bicrystal and Whipple comparison
conda run -n moose ../../../../test/moose_test-opt -i bicrystal_whipple.i
conda run -n moose python scripts/compare_whipple.py

# 3D transient, 1e9 s
conda run -n moose ../../../../test/moose_test-opt -i polycrystal_3d.i
```

To repeat the quadratic 2D experiment on the same coarse topology:

```bash
conda run -n moose ../../../../test/moose_test-opt \
  -i polycrystal_transient.i \
  Mesh/second_order=true \
  Variables/c/order=SECOND \
  AuxVariables/bulk_diffusion_residual/order=SECOND \
  AuxVariables/gb_diffusion_residual/order=SECOND \
  AuxVariables/storage_residual/order=SECOND \
  Outputs/file_base=polycrystal_2d_quadratic
```

This converts the complete mixed-dimensional mesh together, producing shared-node TRI6 and EDGE3
elements. In the tested early transient it eliminated the coarse TRI3 concentration undershoot.

## Regenerate the meshes

The scripts prefer `neper` and `gmsh` on `PATH`, then fall back to the repository-local
`.neper-env/bin` installation. Neper 5.0.0 and Gmsh 4.12.2 were used here.

```bash
# 2D: characteristic length 2.5 um (about four elements per 10 um grain)
./scripts/generate_neper_polycrystal.sh

# Optional 2D refinement
./scripts/generate_neper_polycrystal.sh 1.25 neper_polycrystal_fine

# 3D: characteristic length 5 um
./scripts/generate_neper_polycrystal_3d.sh
```

The exact 3D commands are

```text
neper -T -n 100 -dim 3 -domain 'cube(200,50,50)' -morpho voronoi \
  -id 1 -regularization 1 -o mesh/neper_polycrystal_3d
neper -M mesh/neper_polycrystal_3d.tess -gmsh /absolute/path/to/gmsh \
  -tmp /absolute/path/to/mesh -elttype tet -order 1 -cl 5 \
  -format msh:ascii -o mesh/neper_polycrystal_3d_raw
python3 scripts/preprocess_neper_mesh_3d.py \
  mesh/neper_polycrystal_3d_raw.msh mesh/neper_polycrystal_3d.msh \
  --stats mesh/neper_polycrystal_3d_stats.json
```

Absolute temporary and Gmsh paths matter for 3D Neper meshing. The preprocessors remove points
and edges, combine grain volumes into `bulk`, combine internal faces into `grain_boundaries`, and
label the exterior faces. Per-grain geometrical IDs remain in the Gmsh records.

With 100 grains in the 3D volume, each grain has mean volume 5,000 um^3, or cube-root length
17.1 um. Keeping a 10 um characteristic size in this same 3D domain would require roughly 500
grains.

## Output and checks

CSV records average concentration, bulk inventory, inlet and outlet reaction fluxes, and assembled
mass-balance error. Exodus records both element dimensions for ParaView. Useful finalized files
are `polycrystal_2d.e`, `polycrystal_2d_quadratic.e`, `bicrystal_whipple.e`, and
`polycrystal_3d.e`.

The 2D steady zero-excess case recovers the linear homogeneous solution. Enhanced GB transport
increases its steady total flux by 36.45 times; coarse and refined enhanced fluxes differ by only
0.023%. The Whipple bicrystal relative L2 error falls from `4.74e-4` at 62.5 Ms to `4.40e-5` at
250 Ms.

The 3D output has 46,237 TET10 bulk elements, 9,352 TRI6 GB elements, and 67,199 nodes. At `1e9 s`
its average concentration is `0.3283657` and its outlet flux is `1.579265e-23 m^3/s`. The maximum
assembled balance error is `7.58e-34 m^3/s`. At the first saved state, 100 Ms, eight nodes have a
small Galerkin undershoot with minimum `-2.33e-5`; every saved state from 200 Ms onward is bounded
by zero and one.

## MFEM status

MFEM was intentionally left out of the finalized example. Its dependency builds, but this checkout
unconditionally compiles its `MFEMMUMPS` wrapper while the local conda toolchain cannot link the
Fortran/MUMPS build without the missing macOS SDK. The appropriate future MFEM analogue is one H1
field with `MFEMDiffusionKernel` in the volume and a boundary bilinear `DiffusionIntegrator` on
internal-interface attributes. It should be added only when that backend can be compiled and
tested end to end.
