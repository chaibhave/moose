# Transient diffusion in a 200 um x 50 um x 50 um, 100-grain polycrystal.
length = ${units 200 mum -> m}
center = ${units 25 mum -> m}
D_bulk = ${units 1e-19 m^2/s -> m^2/s}
D_gb = ${fparse 1e6 * D_bulk}
gb_thickness = ${units 0.5 nm -> m}
# A 2D GB embedded in 3D needs K_gb_excess in m^3/s.
K_gb_excess = ${fparse gb_thickness * (D_gb - D_bulk)}
mesh_coordinate_scale = ${units 1 mum -> m}

[Mesh]
  # First-order TET4/TRI3 bulk + TRI3/EDGE2 GB. The early-time undershoot is
  # controlled via mass lumping on the storage kernel (see [Kernels]) instead
  # of promoting the whole mesh to quadratic elements.
  [neper]
    type = FileMeshGenerator
    file = mesh/neper_polycrystal_3d.msh
  []
  [si_units]
    type = TransformGenerator
    input = neper
    transform = SCALE
    vector_value = '${mesh_coordinate_scale} ${mesh_coordinate_scale} ${mesh_coordinate_scale}'
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

[AuxVariables]
  [bulk_diffusion_residual]
    family = LAGRANGE
    order = FIRST
    block = 'bulk grain_boundaries'
  []
  [gb_diffusion_residual]
    family = LAGRANGE
    order = FIRST
    block = 'bulk grain_boundaries'
  []
  [storage_residual]
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
  [gb_transport]
    type = GenericConstantMaterial
    block = grain_boundaries
    prop_names = K_gb_excess_property
    prop_values = ${K_gb_excess}
  []
[]

[Kernels]
  [storage]
    type = MassLumpedTimeDerivative
    variable = c
    block = bulk
    save_in = storage_residual
  []
  [bulk_diffusion]
    type = MatDiffusion
    variable = c
    block = bulk
    diffusivity = D_bulk_property
    save_in = bulk_diffusion_residual
  []
  [gb_excess_diffusion]
    type = MatDiffusion
    variable = c
    block = grain_boundaries
    diffusivity = K_gb_excess_property
    save_in = gb_diffusion_residual
  []
[]

[ICs]
  [initially_empty]
    type = ConstantIC
    variable = c
    value = 0
  []
[]

[BCs]
  [constant_source]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1
  []
  [sink]
    type = DirichletBC
    variable = c
    boundary = right
    value = 0
  []
[]

[Postprocessors]
  [average_c]
    type = ElementAverageValue
    variable = c
    block = bulk
  []
  [bulk_inventory]
    type = ElementIntegralVariablePostprocessor
    variable = c
    block = bulk
  []
  [assembled_storage_rate]
    type = NodalSum
    variable = storage_residual
    block = bulk
  []
  [left_bulk_reaction]
    type = NodalSum
    variable = bulk_diffusion_residual
    boundary = left
    outputs = none
  []
  [left_gb_reaction]
    type = NodalSum
    variable = gb_diffusion_residual
    boundary = left
    outputs = none
  []
  [left_storage_reaction]
    type = NodalSum
    variable = storage_residual
    boundary = left
    outputs = none
  []
  [right_bulk_reaction]
    type = NodalSum
    variable = bulk_diffusion_residual
    boundary = right
    outputs = none
  []
  [right_gb_reaction]
    type = NodalSum
    variable = gb_diffusion_residual
    boundary = right
    outputs = none
  []
  [right_storage_reaction]
    type = NodalSum
    variable = storage_residual
    boundary = right
    outputs = none
  []
  [source_influx]
    type = ParsedPostprocessor
    pp_names = 'left_bulk_reaction left_gb_reaction left_storage_reaction'
    expression = 'left_bulk_reaction+left_gb_reaction+left_storage_reaction'
  []
  [right_outward_flux]
    type = ParsedPostprocessor
    pp_names = 'right_bulk_reaction right_gb_reaction right_storage_reaction'
    expression = '-right_bulk_reaction-right_gb_reaction-right_storage_reaction'
  []
  [conservation_error]
    type = ParsedPostprocessor
    pp_names = 'source_influx right_outward_flux assembled_storage_rate'
    expression = 'source_influx-right_outward_flux-assembled_storage_rate'
  []
[]

[VectorPostprocessors]
  [centerline_profile]
    type = LineValueSampler
    variable = c
    start_point = '0 ${center} ${center}'
    end_point = '${length} ${center} ${center}'
    num_points = 201
    sort_by = id
    execute_on = timestep_end
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  scheme = bdf2
  solve_type = NEWTON
  end_time = ${units 1e9 s -> s}
  dt = ${units 1e7 s -> s}
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  file_base = polycrystal_3d
  csv = true
  exodus = true
  time_step_interval = 10
[]
