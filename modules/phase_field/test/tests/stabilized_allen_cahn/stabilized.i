l = ${units 100.0 m}
V = -500
r0 = ${units 250.0 m}
pi = 3.14159265359
[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = ${fparse 4 * 1000 / ${l} }
    ny = ${fparse 8 * 1000 / ${l} }
    xmin = 0
    xmax = 1000
    ymin = 0
    ymax = 1000
[]
[GlobalParams]
  derivative_order = 2
[]
[Variables]
    [eta]
    []
    [grad_eta]
    []
[]

[Kernels]
    [dt_eta]
        type = ADTimeDerivative
        variable = eta
    []
    [double_well]
        type = ADAllenCahn
        f_name = F
        variable = eta
    []
    [gradient_penalty]
        type = ADACInterface
        variable = eta
        kappa_name = kappa_int
    []
    [stabilize_moelans]
        type = ADACStabilize
        variable = eta
        v = grad_eta
        del_kappa_name = 'del_kappa'
        save_in = 'stab_residual'
        thresh = 1e-8
    []
    [grad_eta_magnitude]
        type = ADGradientMagnitude
        v = eta
        variable = grad_eta
    []
[]
[AuxVariables]
    [stab_residual]
    []
[]

[ICs]
    [eta_IC]
        type = FunctionIC
        function = eta_ic_func
        variable = eta
    []
[]
[Functions]
    [radial_distance]
        type = ParsedFunction
        expression = 'sqrt(x^2+y^2) - ${r0}'
    []
    [eta_ic_func]
        type = ParsedFunction
        symbol_names = 'radial_distance'
        symbol_values = 'radial_distance'
        expression = '0.5-0.5*tanh(2*radial_distance/${l})'
    []
[]

[Materials]
    [Constants]
        type = ADGenericConstantMaterial
        prop_names = 'sigma_int l'
        prop_values = '${units 1 J/m^2} ${units ${l} m}'
    []
    [L]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        property_name = 'L'
        expression = '1'
    []
    [f0]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        property_name = 'f0'
        expression = 'eta^2*(1-eta)^2'
    []
    [kappa_int]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'sigma_int l'
        property_name = 'kappa_int'
        expression = '3*sigma_int*l/4'
    []
    [sigma_s]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'sigma_int l f_chem(eta) mu:=D[f_chem,eta]'
        property_name = 'sigma_s'
        expression = 'max(sigma_int,abs(mu*l/2) )'
        outputs = 'exodus'
    []
    [mu_s]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'sigma_s l'
        property_name = 'mu_s'
        expression = 'sigma_s*6/l'
    []
    [del_kappa]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'sigma_s sigma_int l'
        property_name = 'del_kappa'
        expression = '-3*(sigma_s-sigma_int)*l/4'
        outputs = 'exodus'
    []
    [h]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        property_name = 'h'
        expression = '3*eta^2-2*eta^3'
    []
    [f_chem]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'l h(eta)'
        property_name = f_chem
        expression = 'h*${V}/1.5/l'
        outputs = 'exodus'
    []
    [F]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta grad_eta'
        material_property_names = 'mu_s f0(eta) f_chem(eta)'
        property_name = F
        expression = 'mu_s*f0 + f_chem'
    []
[]

[Postprocessors]
  [total_eta]
    type = ElementIntegralVariablePostprocessor
    variable = eta
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [radius]
    type = ParsedPostprocessor
    expression = 'sqrt(total_eta*4/${pi})'
    pp_names = total_eta
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [expected_radius]
    type = ParsedPostprocessor
    expression = '${r0} - ${V}*t'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [percent_error]
    type = ParsedPostprocessor
    expression = 'abs(expected_radius-radius)*100/expected_radius'
    pp_names = 'expected_radius radius'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Preconditioning]
    [full]
      type = SMP
      full = true
    []
  []
[Executioner]
    type = Transient
    solve_type = PJFNK
    scheme = bdf2
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    # petsc_options_value = 'lu superlu_dist'
    nl_max_its = 10
    # line_search = none
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1e-5
        iteration_window = 2
        optimal_iterations = 7
        growth_factor = 1.25
        cutback_factor = 0.8
      []
    end_time = ${fparse abs( 1.5 * r0 / ${V} ) }
    # num_steps = 1
[]
[Debug]
    show_var_residual_norms = true
[]

[Outputs]
    exodus = true
    csv = true
[]