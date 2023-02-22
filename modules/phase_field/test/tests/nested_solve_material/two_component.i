#
# KKS toy problem in the split form
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 15
  ny = 15
  nz = 0
  xmin = 0
  xmax = 10
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  # elem_type = QUAD4
[]

[Variables]
  # order parameter
  [eta]
    order = FIRST
    family = LAGRANGE
  []
  # hydrogen concentration
  [c_Ni]
  []
  [c_Cr]
  []

[]
[AuxVariables]
  [Fglobal]
    order = CONSTANT
    family = MONOMIAL
  []
    # # chemical potential
  # [w_Ni]
  #   order = FIRST
  #   family = LAGRANGE
  # []

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
    int_width = 0.75
  []
  [c_Ni]
    variable = c_Ni
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = ${fparse 1 - 0.0561 - 2.25398947448865e-07}
    outvalue = 1e-6
    int_width = 0.75
  []
  [c_Cr]
    variable = c_Cr
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 5
    invalue = 0.0561
    outvalue = 1e-4
    int_width = 0.75
  []
[]
# [BCs]
#   [Periodic]
#     [all]
#       variable = 'eta w c'
#       auto_direction = 'x y'
#     []
#   []
# []
[Materials]
  # Free energy of the matrix
  [f_metal]
    type = DerivativeParsedMaterial
    property_name = f_metal
    expression = '96500*(x_cr_m-0.2)^2 + 96500*(x_ni_m-0.8)^2'
    # expression = '8090.0085*x_cr_m*log(x_cr_m) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*log(x_ni_m) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*log(-x_cr_m - x_ni_m + 1) + 123820.210087488'
    # expression = 'x_cr_m*x_ni_m*(10523.5652 - 27907.1252*x_ni_m) + 8090.0085*x_cr_m*log(x_cr_m) - 151569.384260127*x_cr_m + 8090.0085*x_ni_m*log(x_ni_m) - 166834.442737411*x_ni_m + 8090.0085*(-x_cr_m - x_ni_m + 1)*log(-x_cr_m - x_ni_m + 1) + 8090.0085*(5.64178064301999e-79*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^25 + 2.04326786303716e-48*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^15 + 4.89513034938089e-17*(-3605*x_cr_m*x_ni_m - 1109*x_cr_m + 633*x_ni_m)^5)*log(-1.91*x_cr_m*x_ni_m - 2.46*x_cr_m + 0.528*x_ni_m + 1) + 123820.210087488'
    material_property_names = 'x_cr_m x_ni_m'
    additional_derivative_symbols = 'x_cr_m x_ni_m'
    compute = false
  []
  # Free energy of the delta phase
  [f_melt]
    type = DerivativeParsedMaterial
    property_name = f_melt
    expression = '(x_cr_s-1e-4)^2 + (x_ni_s-1e-6)^2'
    # expression = '8089.9721271854*x_cr_s*log(x_cr_s) + x_cr_s*(164050.988348729 - 26181.484*log(973)) + 8089.9721271854*x_ni_s*log(x_ni_s) + x_ni_s*(253092.124999596 - 21499.408*log(973)) + 8089.9721271854*(-x_cr_s - x_ni_s + 1)*log(-x_cr_s - x_ni_s + 1)'
    material_property_names = 'x_cr_s x_ni_s'
    additional_derivative_symbols = 'x_cr_s x_ni_s'
    compute = false
  []

  [equipot_ni]
    type = DerivativeParsedMaterial
    property_name = 'equipot_ni'
    material_property_names = 'f_metal x_ni_m x_cr_m f_melt x_ni_s x_cr_s mu_ni_m:=D[f_metal,x_ni_m] mu_ni_s:=D[f_melt,x_ni_s]'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    expression = mu_ni_m #'(mu_ni_m-mu_ni_s)'
    upstream_materials = 'f_metal f_melt'
    compute = false
  []
  [equipot_cr]
    type = DerivativeParsedMaterial
    material_property_names = 'f_metal x_ni_m x_cr_m mu_cr_m:=D[f_metal,x_cr_m] f_melt x_ni_s x_cr_s  mu_cr_s:=D[f_melt,x_cr_s]'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    property_name = 'equipot_cr'
    upstream_materials = 'f_metal f_melt'
    expression = mu_cr_m #'(mu_cr_m-mu_cr_s)'
    compute = false
  []
  [ni_global_conc]
    type = DerivativeParsedMaterial
    material_property_names = 'x_ni_m x_cr_m x_ni_s x_cr_s h'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr'
    property_name = 'ni_global_conc'
    expression = 'c_Ni - (h*x_ni_m + (1-h)*x_ni_s)'
    compute = false
  []
  [cr_global_conc]
    type = DerivativeParsedMaterial
    material_property_names = 'x_ni_m x_cr_m x_ni_s x_cr_s h'
    additional_derivative_symbols = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    coupled_variables = 'c_Ni c_Cr'
    property_name = 'cr_global_conc'
    expression = 'c_Cr - (h*x_cr_m + (1-h)*x_cr_s)'
    compute = false
  []

  [C]
    type = ParsedMaterial
    material_property_names = 'x_ni_m x_cr_m x_ni_s x_cr_s'
    expression = '(x_ni_m>0)&(x_cr_m>0)&(x_ni_m+x_cr_m<1)' #&(x_ni_s>0)&(x_cr_s>0)
    property_name = 'C'
    compute = false
  []

  # [NestedNewtonSolve]
  #   type = NestedSolveMaterial
  #   xi_names = 'x_ni_m x_cr_m x_ni_s x_cr_s'
  #   Ri = 'ni_global_conc cr_global_conc equipot_ni equipot_cr'
  #   xi_IC = '0.33 0.33 0.33 0.33'
  #   outputs = exodus
  #   output_properties = 'x_ni_m x_cr_m x_ni_s x_cr_s'
  #   absolute_tolerance = 1e-12
  #   relative_tolerance = 1e-8
  #   min_iterations = 1
  #   max_iterations = 5
  # []

  # [TestDampedNewtonSolve]
  #   type = DampedNestedSolveMaterial

  #   xi_names = 'x_ni_m x_cr_m x_ni_s x_cr_s'
  #   Ri = 'ni_global_conc cr_global_conc equipot_ni equipot_cr'
  #   xi_IC = '0.8 0.19 0.33 0.33'
  #   outputs = exodus
  #   output_properties = 'x_ni_m x_cr_m x_ni_s x_cr_s'
  #   absolute_tolerance = 1e-8 #1e-12
  #   relative_tolerance = 1e-8
  #   min_iterations = 1
  #   max_iterations = 50
  #   damping_factor = 0.6
  #   damping_algorithm = BOUNDED_DAMP
  #   conditions = C
  # []
  [TestDampedNewtonSolve]
    type = DampedNestedSolveMaterial

    xi_names = 'x_ni_m x_cr_m'
    Ri = 'equipot_ni equipot_cr'
    xi_IC = '0.8 0.19'
    outputs = exodus
    output_properties = 'x_ni_m x_cr_m'
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-8
    min_iterations = 1
    max_iterations = 50
    damping_factor = 0.6
    damping_algorithm = BOUNDED_DAMP
    conditions = C
    max_damping_iters = 50
  []
  
  # h(eta)
  [h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
    outputs = exodus
    output_properties = 'h'
  []
  # [h_eta]
  #   type = DerivativeParsedMaterial
  #   coupled_variables = 'eta'
  #   expression = 'eta'
  #   property_name = 'h'
  # []
  # g(eta)
  [g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  []
  # constant properties
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M   L   kappa'
    prop_values = '0.7 0.7 0.4'
  []
[]
[Kernels]
  # Cahn-Hilliard Equation
  #
  [dc_Ni_dt]
    type = TimeDerivative
    variable = c_Ni
  []

  [dc_Cr_dt]
    type = TimeDerivative
    variable = c_Cr
  []

  [detadt]
    type = TimeDerivative
    variable = eta
  []
[]
# [AuxKernels]
#   [GlobalFreeEnergy]
#     variable = Fglobal
#     type = KKSGlobalFreeEnergy
#     fa_name = f_metal
#     fb_name = f_melt
#     w = 0.4
#   []
# []
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  # petsc_options_iname = '-pctype -sub_pc_type -sub_pc_factor_shift_type -pc_factor_shift_type'
  # petsc_options_value = ' asm    lu          nonzero                    nonzero'
  l_max_its = 100
  nl_max_its = 100
  num_steps = 1
  dt = 0.1
[]
#
# Precondition using handcoded off-diagonal terms
#
# [Preconditioning]
#   [full]
#     type = SMP
#     full = true
#   []
# []
[Outputs]
  file_base = kks_example_nested
  exodus = true
[]
