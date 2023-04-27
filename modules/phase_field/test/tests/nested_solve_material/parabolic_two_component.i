#
# KKS toy problem in the non-split form
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 80
  ny = 80
  nz = 0
  xmin = -10
  xmax = 10
  ymin = -10
  ymax = 10
  zmin = 0
  zmax = 0
  # elem_type = QUAD4
[]
[Variables]
  # order parameter
  [eta]
    #   order = THIRD
    #   family = HERMITE
  []
  [eta2]

  []
  # hydrogen concentration
  [c]
    #   order = THIRD
    #   family = HERMITE
  []
  [w]

  []

[]
[ICs]
  [eta]
    variable = eta
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 0
    outvalue = 1
    int_width = 2
  []
  [eta2]
    variable = eta2
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 1
    outvalue = 0
    int_width = 2
  []
  [c]
    variable = c
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 0.9 ##0.6
    outvalue = 0.1 ##0.4
    int_width = 2
  []
[]
[BCs]
  # [c_BC]
  #   type = DirichletBC
  #   variable = c
  #   value = 0.8
  #   boundary = 'left right'
  # []
  [w_BC]
    type = DirichletBC
    variable = w
    value = -0.2
    boundary = 'left right'
  []
#   [Periodic]
#     [all]
#       variable = 'eta c'
#       auto_direction = 'x' # y'
#     []
#   []
[]
[Materials]
  # Free energy of the matrix
  [fm]
    type = DerivativeParsedMaterial
    property_name = fm
    material_property_names = 'cm'
    additional_derivative_symbols = 'cm'
    expression = '(0.1-cm)^2'
    compute = false
    outputs = exodus
  []
  # [omega_m]
  #   type = DerivativeParsedMaterial
  #   property_name = omega_m
  #   material_property_names = 'cm fm'
  #   coupled_variables = 'c w eta'
  #   additional_derivative_symbols = 'cm'
  #   expression = 'fm-cm*w'
  #   #outputs = exodus
  #   output_properties = omega_m
  #   compute = false
  # []
  # Free energy of the delta phase
  [fd]
    type = DerivativeParsedMaterial
    property_name = fd
    material_property_names = 'cd'
    additional_derivative_symbols = 'cd'
    expression = '(0.9-cd)^2'
    compute = false
  []
  # [omega_d]
  #   type = DerivativeParsedMaterial
  #   property_name = omega_d
  #   material_property_names = 'cd fd'
  #   coupled_variables = 'c w eta'
  #   additional_derivative_symbols = 'cd'
  #   expression = 'fd-cd*w'
  #   #outputs = exodus
  #   output_properties = omega_d
  #   compute = false
  # []

  # [equipot_c]
  #   type = DerivativeParsedMaterial
  #   property_name = 'equipot_c'
  #   material_property_names = 'fm fd cm cd w_m:=D[fm,cm] w_d:=D[fd,cd]'
  #   additional_derivative_symbols = 'cm cd'
  #   expression = 'w_m - w_d'
  #   upstream_materials = 'fm fd omega_m omega_d'
  #   compute = false
  # []
  # [global_c]
  #   type = DerivativeParsedMaterial
  #   property_name = 'global_c'
  #   material_property_names = 'cm cd h'
  #   coupled_variables = 'c eta'
  #   additional_derivative_symbols = 'cm cd'
  #   expression = 'c - ( h*cm + (1-h)*cd)'
  #   upstream_materials = 'fm fd'
  #   compute = false
  # []

    # [NestedNewtonSolve]
    #   type = NestedSolveMaterial
    #   xi_names = 'cm cd'
    #   Ri = 'equipot_c global_c'
    #   xi_IC = '0.0 0.0'
    #   absolute_tolerance = 1e-12
    #   relative_tolerance = 1e-8
    #   min_iterations = 1
    #   max_iterations = 25
    #   #outputs = exodus
    # []
  # [NestedNewtonSolve]
  #   type = DampedNestedSolveMaterial
  #   xi_names = 'cm cd'
  #   Ri = 'equipot_c global_c'
  #   xi_IC = '0.0 0.0'
  #   absolute_tolerance = 1e-12
  #   relative_tolerance = 1e-8
  #   min_iterations = 1
  #   max_iterations = 25
  #   conditions = C
  #   damping_algorithm = DAMP_FALSE
  #   #outputs = exodus
  # []

  [KKS_nested_solve]
    type = NestedKKSMultiPhaseMaterial
    Fj_material = 'fm fd'
    global_cs = 'c'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'cm cd'
    ci_IC = '0.1 0.1'
    damping_algorithm = BOUNDED_DAMP
    conditions = C
    args = 'c eta'
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-08
    outputs = exodus
  []

  [KKS_nest_conc_derivatives]
    type = KKSPhaseConcentrationMultiPhaseDerivatives
    Fj_material = 'fm fd'
    global_cs = 'c'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'cm cd'
    # outputs = exodus
  []

  [C]
    type = DerivativeParsedMaterial
    property_name = 'C'
    expression = '1.0'
    upstream_materials = 'fd fm'
    compute = false
  []

  # h(eta)
  [h]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = 'h'
    phase_etas = eta
    all_etas = 'eta eta2'
  []

  [h2]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = 'h2'
    phase_etas = eta2
    all_etas = 'eta eta2'
  []
  # g(eta)
  [g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
    #outputs = exodus
  []
  # constant properties
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L M kappa'
    prop_values = '0.7 0.7 0.4' 
  []
[]
[Kernels]
  #
  # Cahn-Hilliard Equation
  #
  [chempot]
    type = NestKKSSplitCHCRes
    global_cs = 'c'
    all_etas = 'eta eta2'
    c1_names = 'cm cd'
    F1_name = 'fm'
    w = w
    variable = c
  []

  [dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
  []
  [MatDif_Ni]
    type = SplitCHWRes
    mob_name = M
    variable = w
    coupled_variables = 'c eta eta2'
  []
  #
  # Allen-Cahn Equation
  #
 
  [detadt]
    type = TimeDerivative
    variable = eta
  []
  [ACBulkF]
    type = NestKKSMultiACBulkF
    variable = eta
    eta_i = eta
    global_cs = 'c'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    Fj_names = 'fm fd'
    ci_names = 'cm cd'
    coupled_variables = 'c w'
    wi = 0.4
    gi_name = g
  []
  [ACBulkC]
    type = NestKKSMultiACBulkC
    variable = eta
    eta_i = eta
    Fj_names = 'fm fd'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'cm cd'
    global_cs = 'c'
    coupled_variables = 'w'
  []
  [ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa
    mob_name = L
  []

  [deta2_dt]
    type = TimeDerivative
    variable = eta2
  []
  [ACBulkF2]
    type = NestKKSMultiACBulkF
    variable = eta2
    eta_i = eta2
    global_cs = 'c'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    Fj_names = 'fm fd'
    ci_names = 'cm cd'
    coupled_variables = 'c w'
    wi = 0.4
    gi_name = g
  []
  [ACBulkC2]
    type = NestKKSMultiACBulkC
    variable = eta2
    eta_i = eta2
    Fj_names = 'fm fd'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'cm cd'
    global_cs = 'c'
    coupled_variables = 'w'
  []
  [ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
    mob_name = L
  []
  
[]

[Postprocessors]
  [total_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Executioner]
  type = Transient
  solve_type = 'NEWTON'
    petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu' 
  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-8
  # num_steps = 100
  end_time = 100
  dt = 0.01
  dtmin = 0.01
[]
[Preconditioning]
  [mydebug]
    type = SMP
    full = true
  []
[]
[Outputs]
  file_base = kks_example
  exodus = true
  csv = true
[]
