# MooseUnits documents each dimensional input and converts it to SI.
length = ${units 200 mum -> m}
D_bulk = ${units 1e-19 m^2/s -> m^2/s}
D_gb = ${fparse 1e6 * D_bulk}
gb_thickness = ${units 0.5 nm -> m}
# K_gb_excess is in m^3/s for a 1D interface embedded in 2D.
K_gb_excess = ${fparse gb_thickness * (D_gb - D_bulk)}
mesh_coordinate_scale = ${units 1 mum -> m}
mesh_file = mesh/neper_polycrystal.msh

[Mesh]
  [neper]
    type = FileMeshGenerator
    file = ${mesh_file}
  []
  [si_units]
    type = TransformGenerator
    input = neper
    transform = SCALE
    vector_value = '${mesh_coordinate_scale} ${mesh_coordinate_scale} 1'
  []
  [explicit_grain_boundaries]
    type = LowerDBlockFromSidesetGenerator
    input = si_units
    sidesets = grain_boundaries
    deduplicate = true
    new_block_id = 1001
    new_block_name = grain_boundaries
  []
[]

[Variables]
  [c]
    family = LAGRANGE
    order = FIRST
    block = 'bulk grain_boundaries'
  []
[]

# These variables collect the assembled bulk and GB residuals. Their sums on
# Dirichlet nodes are the reactions, including point fluxes from GB endpoints.
[AuxVariables]
  [bulk_residual]
    family = LAGRANGE
    order = FIRST
    block = 'bulk grain_boundaries'
  []
  [gb_residual]
    family = LAGRANGE
    order = FIRST
    block = 'bulk grain_boundaries'
  []
[]

[Materials]
  [bulk_transport]
    type = GenericConstantMaterial
    block = bulk
    prop_names = D_bulk_property
    prop_values = ${D_bulk}
  []
  [gb_excess_transport]
    type = GenericConstantMaterial
    block = grain_boundaries
    prop_names = K_gb_excess_property
    prop_values = ${K_gb_excess}
  []
[]

[Kernels]
  [bulk_diffusion]
    type = MatDiffusion
    variable = c
    block = bulk
    diffusivity = D_bulk_property
    save_in = bulk_residual
  []
  # On an embedded EDGE2, libMesh maps the physical shape-function gradient
  # into the edge tangent. MatDiffusion therefore assembles grad_Gamma(c).
  [gb_excess_diffusion]
    type = MatDiffusion
    variable = c
    block = grain_boundaries
    diffusivity = K_gb_excess_property
    save_in = gb_residual
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1
  []
  [right]
    type = DirichletBC
    variable = c
    boundary = right
    value = 0
  []
[]

[Functions]
  [homogeneous_linear_solution]
    type = ParsedFunction
    expression = '1-x/${length}'
  []
[]

[Postprocessors]
  [mesh_elements]
    type = NumElements
    elem_filter = total
  []
  [mesh_nodes]
    type = NumNodes
  []
  [average_c]
    type = ElementAverageValue
    variable = c
    block = bulk
  []
  [linear_profile_l2_error]
    type = ElementL2Error
    variable = c
    function = homogeneous_linear_solution
    block = bulk
  []
  [left_bulk_reaction]
    type = NodalSum
    variable = bulk_residual
    boundary = left
    outputs = none
  []
  [left_gb_reaction]
    type = NodalSum
    variable = gb_residual
    boundary = left
    outputs = none
  []
  [right_bulk_reaction]
    type = NodalSum
    variable = bulk_residual
    boundary = right
    outputs = none
  []
  [right_gb_reaction]
    type = NodalSum
    variable = gb_residual
    boundary = right
    outputs = none
  []
  [left_outward_flux]
    type = ParsedPostprocessor
    pp_names = 'left_bulk_reaction left_gb_reaction'
    expression = '-left_bulk_reaction-left_gb_reaction'
  []
  [right_outward_flux]
    type = ParsedPostprocessor
    pp_names = 'right_bulk_reaction right_gb_reaction'
    expression = '-right_bulk_reaction-right_gb_reaction'
  []
  [flux_imbalance]
    type = ParsedPostprocessor
    pp_names = 'left_outward_flux right_outward_flux'
    expression = 'left_outward_flux+right_outward_flux'
  []
  [right_bulk_flux]
    type = ParsedPostprocessor
    pp_names = right_bulk_reaction
    expression = '-right_bulk_reaction'
  []
  [right_gb_excess_flux]
    type = ParsedPostprocessor
    pp_names = right_gb_reaction
    expression = '-right_gb_reaction'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  automatic_scaling = true
  # A Newton residual evaluation after the linear update populates save_in with
  # the converged reactions used by the boundary-flux postprocessors.
  solve_type = NEWTON
  l_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
  csv = true
[]
