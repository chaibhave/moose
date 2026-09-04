# MooseUnits documents each dimensional input and converts it to SI. The
# bicrystal is long in x, with one straight GB at y = width / 2.
length = ${units 200 mum -> m}
width = ${units 50 mum -> m}
gb_y = ${units 25 mum -> m}
D_bulk = ${units 1e-19 m^2/s -> m^2/s}
D_gb = ${fparse 1e6 * D_bulk}
gb_thickness = ${units 0.5 nm -> m}
# K_gb_excess is in m^3/s for a 1D interface embedded in 2D.
K_gb_excess = ${fparse gb_thickness * (D_gb - D_bulk)}

[Mesh]
  [rectangle]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 50
    xmin = 0
    xmax = ${length}
    ymin = 0
    ymax = ${width}
    elem_type = TRI3
  []
  [lower_grain]
    type = SubdomainBoundingBoxGenerator
    input = rectangle
    bottom_left = '0 0 0'
    top_right = '${length} ${gb_y} 0'
    block_id = 1
    block_name = lower_grain
  []
  [upper_grain]
    type = SubdomainBoundingBoxGenerator
    input = lower_grain
    bottom_left = '0 ${gb_y} 0'
    top_right = '${length} ${width} 0'
    block_id = 2
    block_name = upper_grain
  []
  [gb_sideset]
    type = SideSetsBetweenSubdomainsGenerator
    input = upper_grain
    primary_block = lower_grain
    paired_block = upper_grain
    new_boundary = grain_boundaries
  []
  [explicit_grain_boundary]
    type = LowerDBlockFromSidesetGenerator
    input = gb_sideset
    sidesets = grain_boundaries
    new_block_id = 3
    new_block_name = grain_boundaries
  []
[]

[Variables]
  [c]
    family = LAGRANGE
    order = FIRST
    block = 'lower_grain upper_grain grain_boundaries'
  []
[]

[AuxVariables]
  [bulk_diffusion_residual]
    family = LAGRANGE
    order = FIRST
    block = 'lower_grain upper_grain grain_boundaries'
  []
  [gb_diffusion_residual]
    family = LAGRANGE
    order = FIRST
    block = 'lower_grain upper_grain grain_boundaries'
  []
  [storage_residual]
    family = LAGRANGE
    order = FIRST
    block = 'lower_grain upper_grain grain_boundaries'
  []
[]

[Materials]
  [bulk_transport]
    type = GenericConstantMaterial
    block = 'lower_grain upper_grain'
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
    type = TimeDerivative
    variable = c
    block = 'lower_grain upper_grain'
    save_in = storage_residual
  []
  [bulk_diffusion]
    type = MatDiffusion
    variable = c
    block = 'lower_grain upper_grain'
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
[]

[Postprocessors]
  [bulk_inventory]
    type = ElementIntegralVariablePostprocessor
    variable = c
    block = 'lower_grain upper_grain'
    execute_on = 'initial timestep_end'
  []
  [inventory_finite_difference]
    type = ChangeOverTimePostprocessor
    postprocessor = bulk_inventory
    divide_by_dt = true
    execute_on = 'initial timestep_end'
  []
  [assembled_storage_residual]
    type = NodalSum
    variable = storage_residual
    block = 'lower_grain upper_grain'
    outputs = none
  []
  [assembled_storage_rate]
    type = ParsedPostprocessor
    pp_names = assembled_storage_residual
    expression = 'assembled_storage_residual'
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
  [source_influx]
    type = ParsedPostprocessor
    pp_names = 'left_bulk_reaction left_gb_reaction left_storage_reaction'
    expression = 'left_bulk_reaction+left_gb_reaction+left_storage_reaction'
  []
  [conservation_error]
    type = ParsedPostprocessor
    pp_names = 'source_influx assembled_storage_rate'
    expression = 'source_influx-assembled_storage_rate'
  []
[]

[VectorPostprocessors]
  [gb_profile]
    type = LineValueSampler
    variable = c
    start_point = '0 ${gb_y} 0'
    end_point = '${units 80 mum -> m} ${gb_y} 0'
    num_points = 161
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
  start_time = 0
  end_time = ${units 2.5e8 s -> s}
  dt = ${units 2.5e6 s -> s}
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-10
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  file_base = bicrystal_whipple
  [csv]
    type = CSV
    time_step_interval = 25
    time_data = true
    precision = 14
  []
  [exodus]
    type = Exodus
    time_step_interval = 10
  []
[]
