#
# KKS toy problem in the split form
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100 #${fparse int ( 10 * ( ${xmax} - ${xmin} ) ) }
  ny = 15
  nz = 0
  xmin = 0
  xmax = 10 #0
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  # elem_type = QUAD4
[]

[GlobalParams]
  # enable_jit = false
  # profile = TANH
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
    order = FIRST
    family = LAGRANGE
    # initial_condition = 0.268941421370010
  []
[]
[AuxVariables]
  [Fglobal]
    order = CONSTANT
    family = MONOMIAL
  []

[]
[ICs]
  [eta]
    variable = eta
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 1
    outvalue = 0
    int_width = 1
  []
  [c_Ni]
    variable = c_Ni
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 0.9499997060782867 ##'${fparse 1 - 0.0561 - 2.25398947448865e-07}'
    outvalue = 0.0000000080108391 ##1e-10
    int_width = 1
  []
  [c_Cr]
    variable = c_Cr
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 0.0500000000000000 ##0.0561
    outvalue = 0.0089901478135457 ##1e-4
    int_width = 1
  []
  # [c_Ni]
  #   variable = c_Ni
  #   type = SmoothCircleIC
  #   x1 = 0.0
  #   y1 = 0.0
  #   radius = 5
  #   invalue = '${fparse 1 - 0.0561 - 2.25398947448865e-07}'
  #   outvalue = 1e-10
  #   int_width = 1
  # []
  # [c_Cr]
  #   variable = c_Cr
  #   type = SmoothCircleIC
  #   x1 = 0.0
  #   y1 = 0.0
  #   radius = 5
  #   invalue = 0.0561
  #   outvalue = 1e-4
  #   int_width = 1
  # []
[]
# [BCs]
#   [c_Ni_right]
#     type = DirichletBC
#     variable = c_Ni
#     value = 8e-9 ##${fparse 1 - 0.0561 - 2.25398947448865e-07} ##
#     boundary = right
#   []
#   [c_Cr_right]
#     type = DirichletBC
#     variable = c_Cr
#     value = 4e-3 ##0.0561 ##
#     boundary = right
#   []
# []
[Materials]
  # Free energy of the matrix
  [f_metal]
    type = DerivativeParsedMaterial
    property_name = f_metal
    expression = '(8090.0085*x_cr_m*log(x_cr_m) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*log(x_ni_m) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*log(-x_cr_m - x_ni_m + 1) + 123820.210087488)/96500'
    # expression = 'x_cr_m*x_ni_m*(10523.5652 - 27907.1252*x_ni_m) + 8090.0085*x_cr_m*log(x_cr_m) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*log(x_ni_m) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*log(-x_cr_m - x_ni_m + 1) + 8090.0085*(5.64178064301999e-79*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^25 + 2.04326786303716e-48*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^15 + 4.89513034938089e-17*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^5)*log(-1.91*x_cr_m*x_ni_m - 2.46*x_cr_m + 0.528*x_ni_m + 1) + 123820.210087488'
    material_property_names = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    additional_derivative_symbols = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    derivative_order = 3
    compute = false
    outputs = exodus
  []
  [omega_metal]
    type = DerivativeParsedMaterial
    property_name = omega_metal
    coupled_variables = 'eta c_Ni c_Cr w_Cr w_Ni'
    material_property_names = 'x_ni_m x_cr_m f_metal'
    expression = 'f_metal - x_ni_m*w_Ni - x_cr_m*w_Cr'
    compute = false
    outputs = exodus
    output_properties = 'omega_metal'
  []
  # Free energy of the delta phase
  [f_melt]
    type = DerivativeParsedMaterial
    property_name = f_melt
    expression = '(8089.9721271854*x_cr_s*log(x_cr_s) + x_cr_s*(164050.988348729 - 26181.484*log(973)) + 8089.9721271854*x_ni_s*log(x_ni_s) + x_ni_s*(253092.124999596 - 21499.408*log(973)) + 8089.9721271854*(-x_cr_s - x_ni_s + 1)*log(-x_cr_s - x_ni_s + 1))/96500'
    material_property_names = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    additional_derivative_symbols = 'x_cr_m x_ni_m x_cr_s x_ni_s'
    derivative_order = 3
    compute = false
    outputs = exodus
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    property_name = omega_melt
    coupled_variables = 'eta c_Ni c_Cr w_Cr w_Ni'
    material_property_names = 'x_ni_s x_cr_s f_melt'
    expression = 'f_melt - x_ni_s*w_Ni - x_cr_s*w_Cr'
    compute = false
    outputs = exodus
    output_properties = 'omega_melt'
  []

  [equipot_ni]
    type = DerivativeParsedMaterial
    property_name = 'equipot_ni'
    material_property_names = 'f_metal(x_ni_m,x_cr_m) x_ni_m x_cr_m mu_ni_m:=D[f_metal,x_ni_m] f_melt(x_ni_s,x_cr_s) x_ni_s x_cr_s mu_ni_s:=D[f_melt,x_ni_s]'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s f_metal f_melt'
    # expression = 'mu_ni_m - mu_ni_s'
    expression = '(8090.0085*log(x_ni_m) - 8089.9721271854*log(x_ni_s) - 8090.0085*log(-x_cr_m - x_ni_m + 1) + 8089.9721271854*log(-x_cr_s - x_ni_s + 1) - 419926.567737007 + 21499.408*log(973))/96500'
    derivative_order = 2
    upstream_materials = 'f_metal f_melt omega_metal omega_melt'
    compute = false
  []
  [equipot_cr]
    type = DerivativeParsedMaterial
    material_property_names = 'f_metal(x_ni_m,x_cr_m) x_ni_m x_cr_m mu_cr_m:=D[f_metal,x_cr_m] f_melt(x_ni_s,x_cr_s) x_ni_s x_cr_s  mu_cr_s:=D[f_melt,x_cr_s]' #mu_cr_m:=D[f_metal,x_cr_m] 
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s f_metal f_melt'
    property_name = 'equipot_cr'
    upstream_materials = 'f_metal f_melt omega_metal omega_melt'
    # expression = 'mu_cr_m - mu_cr_s' 
    expression = '(8090.0085*log(x_cr_m) - 8089.9721271854*log(x_cr_s) - 8090.0085*log(-x_cr_m - x_ni_m + 1) + 8089.9721271854*log(-x_cr_s - x_ni_s + 1) - 315620.372608856 + 26181.484*log(973))/96500'
    derivative_order = 2
    compute = false
  []
  [ni_global_conc]
    type = DerivativeParsedMaterial
    material_property_names = 'x_ni_m x_ni_s x_cr_m x_cr_s h'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr'
    property_name = 'ni_global_conc'
    expression = 'c_Ni - (h*x_ni_m + (1-h)*x_ni_s)'
    derivative_order = 2
    compute = false
  []
  [cr_global_conc]
    type = DerivativeParsedMaterial
    material_property_names = 'x_ni_m x_ni_s x_cr_m x_cr_s h'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr'
    property_name = 'cr_global_conc'
    expression = 'c_Cr - (h*x_cr_m + (1-h)*x_cr_s)'
    derivative_order = 2
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
  []

  [TestDampedNewtonSolve]
    type = DampedNestedSolveMaterial

    xi_names = 'x_ni_m x_ni_s x_cr_m x_cr_s'
    Ri = 'equipot_ni ni_global_conc equipot_cr cr_global_conc'
    xi_IC = '0.33 0.33 0.33 0.33'
    outputs = exodus
    output_properties = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-8
    min_iterations = 1
    max_iterations = 100
    damping_factor = 0.6
    damping_algorithm = BOUNDED_DAMP
    conditions = C
    max_damping_iters = 30
    delta_X_threshold = 2e-16
  []
  
  [dx_ni_m_dc_Ni]
    type = DerivativeParsedMaterial
    coupled_variables = 'eta c_Cr c_Ni'
    material_property_names = 'x_ni_m x_ni_s x_cr_m x_cr_s h dh_deta:=D[h,eta] f_metal f_melt f2_metal:=D[f_metal,x_ni_m,x_ni_m] f2_melt:=D[f_melt,x_ni_s,x_ni_s]'
    expression = 'f2_melt/(h*f2_melt + (1-h)*f2_metal )'
    property_name = 'dx_ni_m/dc_Ni'
  []
  
  # h(eta)
  [h]
    type = DerivativeParsedMaterial
    coupled_variables = 'eta'
    expression = 'eta^2/(eta^2 + (1-eta)^2 )'
    property_name = 'h'
  []
  # g(eta)
  [g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  []
  # constant properties
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M     L   kappa'
    prop_values = '1e-3 1e-4 0.4'
  []
[]

[Kernels]
  # Cahn-Hilliard Equation
  #
  [chempot_Ni]
    type = CHSplitChemicalPotential
    variable = w_Ni
    chemical_potential_prop = 'df_melt/dx_ni_s'
    c = c_Ni
  []

  [dc_Ni_dt]
    type = TimeDerivative
    variable = c_Ni
  []
  [MatDif_Ni]
    type = MatDiffusion
    diffusivity = M
    variable = c_Ni
    v = w_Ni
    args = 'eta w_Cr'
  []

  [chempot_Cr]
    type = CHSplitChemicalPotential
    variable = w_Cr
    chemical_potential_prop = 'df_melt/dx_cr_s'
    c = c_Cr
  []
  [dc_Cr_dt]
    type = TimeDerivative
    variable = c_Cr
  []
  [MatDif_Cr]
    type = MatDiffusion
    diffusivity = M
    variable = c_Cr
    v = w_Cr
    args = 'eta w_Ni'
  []

  [detadt]
    type = TimeDerivative
    variable = eta
  []
  # [ACBulkF]
  #   type = KKSACBulkF
  #   variable = eta
  #   fa_name = omega_metal
  #   fb_name = omega_melt
  #   coupled_variables = 'c_Cr c_Ni w_Cr w_Ni'
  #   w = 0.4
  # []
  # [ACBulkC]
  #   type = KKSACBulkC
  #   variable = eta
  #   ca = c_Ni
  #   cb = cd
  #   fa_name = fm
  #   mob_name = L
  # []
  # [ACInterface]
  #   type = ACInterface
  #   variable = eta
  #   kappa_name = kappa
  #   mob_name = L
  # []
[]

# [UserObjects]
#   [NP_UO_w_Ni]
#     type = NodalPatchRecoveryMaterialProperty
#     patch_polynomial_order = SECOND
#     property = 'df_melt/dx_ni_s'

#     execute_on = 'TIMESTEP_END'
#   []
# []

[AuxKernels]
  #   [GlobalFreeEnergy]
  #     variable = Fglobal
  #     type = KKSGlobalFreeEnergy
  #     fa_name = f_metal
  #     fb_name = f_melt
  #     w = 0.4
  #   []
  # [w_Ni_val]
  #   type = NodalPatchRecoveryAux
  # []

[]
[Executioner]
  type = Transient
  solve_type = NEWTON ##PJFNK
  scheme = bdf2
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_tol = 1e-3
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-9
  l_max_its = 100
  nl_max_its = 100
  end_time = 10
  num_steps = 1 #1000
  # automatic_scaling = true
  # dt = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-6
    iteration_window = 2
    optimal_iterations = 9
    growth_factor = 1.25
    cutback_factor = 0.8
  []
[]

[Debug]
  show_var_residual_norms = true
[]

#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  # [full]
  #   type = SMP
  #   full = true
  # []
  [FD]
    type = FDP
    full = true
  []
[]
[Outputs]
  file_base = nested_solve
  exodus = true
[]
