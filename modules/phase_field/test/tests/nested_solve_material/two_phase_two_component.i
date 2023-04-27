#
# KKS toy problem in the split form
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${fparse int ( 10 * xmax / sqrt(2) ) }
  ny = 15
  nz = 0
  xmin = 0
  xmax = 50 #0
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  # uniform_refine = 1
  # skip_partitioning = true
  # elem_type = QUAD4
[]

[GlobalParams]
  enable_jit = false
  # profile = TANH
  x_Cr_BC = 1e-4
  x_Ni_BC = 8e-9 ##0.0000000080108391

  # mu_Ni_BC = ${fparse ( 8089.9721271854 * log( ${GlobalParams/x_Ni_BC} ) - 8089.9721271854 * log( - ${GlobalParams/x_Cr_BC} - ${GlobalParams/x_Ni_BC} + 1 ) - 21499.408*log(973) + 253092.124999596 ) / 96500 }
  # mu_Cr_BC = ${fparse ( 8089.9721271854*log( ${GlobalParams/x_Cr_BC} ) - 8089.9721271854*log(- ${GlobalParams/x_Cr_BC} - ${GlobalParams/x_Ni_BC} + 1 ) - 26181.484*log(973) + 164050.988348729 ) / 96500 }
[]

[Variables]
  # hydrogen concentration
  [c_Ni]
    # initial_condition = 0.253854478070618
  []
  [c_Cr]
    # initial_condition = 0.0151607195967205
  []
  [w_Ni]
  []
  [w_Cr]
  []
  # order parameter
  [eta]
  []
  [eta2]
  []
[]
[AuxVariables]
  [Fglobal]
    order = CONSTANT
    family = MONOMIAL
  []
  [bnds_dummy]
  []
[]
[AuxKernels]
  # [Fglobal_aux]
  #   type = KKSGlobalFreeEnergy
  #   variable = Fglobal
  #   interfacial_vars = 'eta eta2'
  #   kappa_names = 'kappa kappa'
  #   w = 0.4
  #   fa_name = f_metal
  #   fb_name = f_melt
  #   h_name = 'h h2'
  #   g_name = g
  # []
  [Fglobal_aux]
    type = KKSMultiFreeEnergy
    variable = Fglobal
    interfacial_vars = 'eta eta2'
    kappa_names = 'kappa kappa'
    w = 0.4
    Fj_names = 'f_metal f_melt'
    gj_names = 'g g'
    hj_names = 'h h2'
  []
[]

[ICs]
  [eta]
    variable = eta
    type = FunctionIC
    function = eta_func
  []
  [eta2]
    variable = eta2
    type = FunctionIC
    function = eta2_func
  []
  [c_Ni]
    variable = c_Ni
    type = FunctionIC
    function = c_Ni_func
  []
  [c_Cr]
    variable = c_Cr
    type = FunctionIC
    function = c_Cr_func
  []
[]

[Functions]
  [eta_func]
    type = ParsedFunction
    expression = '0.5*(1-tanh((x- ${fparse ${Mesh/xmax} / 2 } )/sqrt(2) ))'
  []
  [eta2_func]
    type = ParsedFunction
    expression = '0.5*(1+tanh((x- ${fparse ${Mesh/xmax} / 2 } )/sqrt(2) ))'
  []
  [h_func]
    type = ParsedFunction
    symbol_names = 'N1 N2'
    symbol_values = 'eta_func eta2_func'
    expression = 'N1^2/(N1^2 + N2^2)'
  []
  [c_Ni_func]
    type = ParsedFunction
    symbol_names = 'h1 '
    symbol_values = 'h_func'
    expression = 'h1*0.9499997060782867 + (1-h1)*0.0000000080108391'
  []
  [c_Cr_func]
    type = ParsedFunction
    symbol_names = 'h1 '
    symbol_values = 'h_func'
    expression = 'h1*0.05 + (1-h1)*0.0089901478135457'
  []

  [c_Ni_BC_func]
    type = ParsedFunction
    expression = '8.0e-9 - (8.0e-9 - ${GlobalParams/x_Ni_BC} )'
  []
  [c_Cr_BC_func]
    type = ParsedFunction
    expression = '0.0000000080108391 - (0.0000000080108391 - ${GlobalParams/x_Cr_BC})'
  []
  [w_Ni_BC_func]
    type = ParsedFunction
    symbol_names = 'x_ni_s x_cr_s'
    symbol_values = 'c_Ni_BC_func c_Cr_BC_func'
    expression = '(8089.9721271854*log(x_ni_s) - 8089.9721271854*log(-x_cr_s - x_ni_s + 1) - 21499.408*log(973) + 253092.124999596)/96500'
  []
  [w_Cr_BC_func]
    type = ParsedFunction
    symbol_names = 'x_ni_s x_cr_s'
    symbol_values = 'c_Ni_BC_func c_Cr_BC_func'
    expression = '(8089.9721271854*log(x_cr_s) - 8089.9721271854*log(-x_cr_s - x_ni_s + 1) - 26181.484*log(973) + 164050.988348729)/96500'
  []

[]

[BCs]
  [c_Ni_right]
    type = FunctionDirichletBC
    variable = c_Ni
    function = c_Ni_BC_func
    # value = ${GlobalParams/x_Ni_BC} ##${fparse 1 - 0.0561 - 2.25398947448865e-07} ##
    boundary = right
  []
  [c_Cr_right]
    type = FunctionDirichletBC
    variable = c_Cr
    function = c_Cr_BC_func
    # value = ${GlobalParams/x_Cr_BC} ##0.0561 ##
    boundary = right
  []

  [w_Ni_right]
    type = FunctionDirichletBC
    variable = w_Ni
    function = w_Ni_BC_func
    # value = '${GlobalParams/mu_Ni_BC}'
    boundary = right
  []
  [w_Cr]
    type = FunctionDirichletBC
    variable = w_Cr
    function = w_Cr_BC_func
    # value = '${GlobalParams/mu_Cr_BC}'
    boundary = right
  []
[]
[Materials]
  # Free energy of the matrix
  [f_metal]
    type = DerivativeParsedMaterial
    property_name = f_metal
    # expression = '(8090.0085*x_cr_m*log(x_cr_m) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*log(x_ni_m) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*log(-x_cr_m - x_ni_m + 1) + 123820.210087488)/96500'
    expression = '(8090.0085*x_cr_m*plog(x_cr_m,err) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*plog(x_ni_m,err) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*plog(-x_cr_m - x_ni_m + 1,err) + 123820.210087488)/96500'
    # expression = '(x_cr_m*x_ni_m*(10523.5652 - 27907.1252*x_ni_m) + 8090.0085*x_cr_m*plog(x_cr_m,err) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*plog(x_ni_m,err) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*plog(-x_cr_m - x_ni_m + 1,err) + 8090.0085*(5.64178064301999e-79*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^25 + 2.04326786303716e-48*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^15 + 4.89513034938089e-17*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^5)*plog(-1.91*x_cr_m*x_ni_m - 2.46*x_cr_m + 0.528*x_ni_m + 1,err) + 123820.210087488)/96500'
    material_property_names = 'x_cr_m x_ni_m x_cr_s x_ni_s err'
    additional_derivative_symbols = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    derivative_order = 3
    compute = false
    # outputs = exodus
  []

  # Free energy of the delta phase
  [f_melt]
    type = DerivativeParsedMaterial
    property_name = f_melt
    # expression = '(8089.9721271854*x_cr_s*log(x_cr_s) + x_cr_s*(164050.988348729 - 26181.484*log(973)) + 8089.9721271854*x_ni_s*log(x_ni_s) + x_ni_s*(253092.124999596 - 21499.408*log(973)) + 8089.9721271854*(-x_cr_s - x_ni_s + 1)*log(-x_cr_s - x_ni_s + 1))/96500'
    expression = '(8089.9721271854*x_cr_s*plog(x_cr_s,err) + x_cr_s*(164050.988348729 - 26181.484*log(973)) + 8089.9721271854*x_ni_s*plog(x_ni_s,err) + x_ni_s*(253092.124999596 - 21499.408*log(973)) + 8089.9721271854*(-x_cr_s - x_ni_s + 1)*plog(-x_cr_s - x_ni_s + 1,err))/96500'
    material_property_names = 'x_cr_m x_ni_m x_cr_s x_ni_s err'
    additional_derivative_symbols = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    derivative_order = 3
    compute = false
    ##outputs = exodus
  []

  [omega_chem]
    type = DerivativeParsedMaterial
    property_name = omega_chem
    material_property_names = 'f_metal f_melt x_ni_m x_cr_m x_ni_s x_cr_s h h2'
    coupled_variables = 'w_Ni w_Cr'
    expression = '(f_metal-w_Ni*x_ni_m-w_Cr*x_cr_m)*h + (f_melt-w_Ni*x_ni_s-w_Cr*x_cr_s)*h2'
    # outputs = exodus
    output_properties = omega_chem
    compute = false
    derivative_order = 0
  []
  [M_Ni]
    type = DerivativeParsedMaterial
    property_name = M_Ni
    material_property_names = 'f_metal f_melt x_ni_m x_cr_m x_ni_s x_cr_s h h2 d2f_metal:=D[f_metal,x_ni_m,x_ni_m] d2f_melt:=D[f_melt,x_ni_s,x_ni_s]'
    expression = '2e-7*h/d2f_metal + 1.0*h2/d2f_melt'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr w_Ni w_Cr eta eta2'
    compute = false
  []
  [M_Cr]
    type = DerivativeParsedMaterial
    property_name = M_Cr
    material_property_names = 'f_metal f_melt x_ni_m x_cr_m x_ni_s x_cr_s h h2 d2f_metal:=D[f_metal,x_cr_m,x_cr_m] d2f_melt:=D[f_melt,x_cr_s,x_cr_s]'
    expression = '2e-7*h/d2f_metal + 1.0*h2/d2f_melt'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr w_Ni w_Cr eta eta2'
    compute = false
  []

  [C]
    type = DerivativeParsedMaterial
    material_property_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    expression = '(x_ni_m>0.0) & (x_ni_s>0.0) & (x_cr_m>0.0) & (x_cr_s>0.0) & ((x_ni_m+x_cr_m)<1.0) & (x_ni_s+x_cr_s<1.0)'
    property_name = 'C'
    derivative_order = 0
    enable_jit = true
    epsilon = 0.0
    disable_fpoptimizer = true
    compute = false
    upstream_materials = 'omega_chem M_Ni M_Cr'
  []

  [KKS_nested_solve]
    type = NestedKKSMultiPhaseMaterial
    Fj_material = 'f_metal f_melt'
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    ci_IC = '0.9 1e-6 0.09 1e-4'
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-8
    min_iterations = 1
    max_iterations = 100
    damping_factor = 0.8
    damping_algorithm = BOUNDED_DAMP
    conditions = C
    max_damping_iters = 50
    delta_X_threshold = 1e-20 #2e-16
    # outputs = exodus
  []

  [KKS_phase_conc_derivatives]
    type = KKSPhaseConcentrationMultiPhaseDerivatives
    Fj_material = 'f_metal f_melt'
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    # outputs = exodus
  []

  [eta_consv_check]
    type = ParsedMaterial
    coupled_variables = 'eta eta2'
    property_name = 'sum_eta'
    expression = '1-eta-eta2'
    # outputs =  exodus
  []

  # h(eta)
  [h]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = 'h'
    phase_etas = eta
    all_etas = 'eta eta2'
    # outputs = exodus
  []

  [h2]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = 'h2'
    phase_etas = eta2
    all_etas = 'eta eta2'
    # outputs = exodus
  []
  # g(eta)
  [g_eta]
    type = MultiBarrierFunctionMaterial ##BarrierFunctionMaterial
    etas = 'eta eta2'
  []
  [g2_eta]
    type = MultiBarrierFunctionMaterial ##BarrierFunctionMaterial
    etas = 'eta eta2'
    # type = BarrierFunctionMaterial
    # g_order = SIMPLE
    # eta = eta2
    function_name = g2
  []
  # constant properties
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M     L   kappa err'
    prop_values = '1e-4 1e-4 0.4 1e-15'
  []
  
[]

[Kernels]
  # Cahn-Hilliard Equation
  #
  [chempot_Ni]
    type = NestKKSSplitCHCRes
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    c1_names = 'x_ni_m x_cr_m'
    F1_name = 'f_metal'
    w = w_Ni
    variable = c_Ni
  []

  [dc_Ni_dt]
    type = CoupledTimeDerivative
    variable = w_Ni
    v = c_Ni
  []
  [MatDif_Ni]
    type = SplitCHWRes
    mob_name = M_Ni
    variable = w_Ni
    coupled_variables = 'c_Ni c_Cr eta eta2 w_Cr'
  []

  [chempot_Cr]
    type = NestKKSSplitCHCRes
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    c1_names = 'x_ni_m x_cr_m'
    F1_name = 'f_metal'
    w = w_Cr
    variable = c_Cr
  []
  [dc_Cr_dt]
    type = CoupledTimeDerivative
    variable = w_Cr
    v = c_Cr
  []
  [MatDif_Cr]
    type = SplitCHWRes
    mob_name = M_Cr
    variable = w_Cr
    coupled_variables = 'c_Ni c_Cr eta eta2 w_Ni'
  []

  [detadt]
    type = TimeDerivative
    variable = eta
  []
  [ACBulkF]
    type = NestKKSMultiACBulkF
    variable = eta
    eta_i = eta
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    Fj_names = 'f_metal f_melt'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    coupled_variables = 'c_Cr c_Ni w_Cr w_Ni'
    wi = 0.4
    gi_name = g
  []
  [ACBulkC]
    type = NestKKSMultiACBulkC
    variable = eta
    eta_i = eta
    Fj_names = 'f_metal f_melt'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    global_cs = 'c_Ni c_Cr'
    coupled_variables = 'w_Ni w_Cr'
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
    global_cs = 'c_Ni c_Cr'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    Fj_names = 'f_metal f_melt'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    coupled_variables = 'c_Cr c_Ni w_Cr w_Ni'
    wi = 0.4
    gi_name = g2
  []
  [ACBulkC2]
    type = NestKKSMultiACBulkC
    variable = eta2
    eta_i = eta2
    Fj_names = 'f_metal f_melt'
    all_etas = 'eta eta2'
    hj_names = 'h h2'
    ci_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    global_cs = 'c_Ni c_Cr'
    coupled_variables = 'w_Ni w_Cr'
  []
  [ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
    mob_name = L
  []
[]
# [Bounds]
#   [c_Ni_bnds_upper]
#     type = ConstantBoundsAux
#     variable = bnds_dummy
#     bounded_variable = c_Ni
#     bound_type = upper
#     bound_value = 1.0
#   []
#   [c_Ni_bnds_lower]
#     type = ConstantBoundsAux
#     variable = bnds_dummy
#     bounded_variable = c_Ni
#     bound_type = lower
#     bound_value = 1e-12
#   []
#   [c_Cr_bnds_upper]
#     type = ConstantBoundsAux
#     variable = bnds_dummy
#     bounded_variable = c_Cr
#     bound_type = upper
#     bound_value = 1.0
#   []
#   [c_Cr_bnds_lower]
#     type = ConstantBoundsAux
#     variable = bnds_dummy
#     bounded_variable = c_Cr
#     bound_type = lower
#     bound_value = 1e-12
#   []
# []

[Dampers]
  [c_Ni]
    type = BoundingValueElementDamper
    variable = c_Ni
    max_value = '${fparse 1.0 - 1e-30}'
    min_value = 1e-16
    min_damping = 0.1
  []
  [c_Cr]
    type = BoundingValueElementDamper
    variable = c_Cr
    max_value = '${fparse 1.0 - 1e-30}'
    min_value = 1e-16
    min_damping = 0.1
  []
  # [eta]
  #   type = BoundingValueElementDamper
  #   variable = eta
  #   max_value = '${fparse 1.0 - 1e-30}'
  #   min_value = 1e-16
  #   min_damping = 0.1
  # []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  # petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm       lu            nonzero'

  # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           2'
  # petsc_options_iname = '-pc_type -snes_type'
  # petsc_options_value = 'lu vinewtonssls' 
  # petsc_options_iname = '-pc_type -snes_type -pc_factor_shift_type -pc_factor_shift_amount '
  # petsc_options_value = 'lu vinewtonrsls NONZERO 1e-10'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type' # -pc_factor_shift_type -pc_factor_shift_amount 
  petsc_options_value = 'lu superlu_dist vinewtonrsls' # NONZERO 1e-10
# petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -snes_type'# -pc_factor_shift_type -pc_factor_shift_amount '
#   petsc_options_value = 'hypre    boomeramg      31                 0.7 vinewtonrsls'# NONZERO 1e-10'

  # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap -snes_type'# -pc_factor_shift_type -pc_factor_shift_amount '
  # petsc_options_value = 'asm      31                  preonly       lu           2 vinewtonrsls'# NONZERO 1e-10'

  
  l_tol = 1e-3
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
  l_max_its = 100
  nl_max_its = 18
  end_time = 3.6e6
  # num_steps = 1000
  # automatic_scaling = true
  dt = 1e-3
  dtmin = 1e-20
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    iteration_window = 2
    optimal_iterations = 15
    growth_factor = 1.25
    cutback_factor = 0.8
  []
  # [Adaptivity]
  #   max_h_level = 2
  #   refine_fraction = 0.2
  #   coarsen_fraction = 0.1
  # []
[]
[Postprocessors]
  [elapsed]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
  [total_eta]
    type = ElementIntegralVariablePostprocessor
    variable = eta
    execute_on = 'INITIAL TIMESTEP_END'

  []
  [total_Ni]
    type = ElementIntegralVariablePostprocessor
    variable = c_Ni
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_Cr]
    type = ElementIntegralVariablePostprocessor
    variable = c_Cr
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = Fglobal
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# [Debug]
#   show_var_residual_norms = true
# []

#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [full]
    type = SMP
    full = true
  []
  # [FD]
  #   type = FDP
  #   full = true
  # []
[]
[Outputs]
  file_base = kks_example_nested
  exodus = true
  csv = true
  perf_graph = true
  # minimum_time_interval = 0.1
[]
