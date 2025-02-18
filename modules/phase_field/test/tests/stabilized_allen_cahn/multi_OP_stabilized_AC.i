l = ${units 20.0 m}
V = -5
r0 = ${units 250.0 m}
f_stab = 2
pi = 3.14159265359
x0 = 500
y0 = 700
[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = ${fparse 4 * 1000 / ${l} }
    ny = ${fparse 4 * 1000 / ${l} }
    xmin = 0
    xmax = 1000
    ymin = 0
    ymax = 1000
[]
[GlobalParams]
    derivative_order = 2
  []
[Variables]
    [eta1]
    []
    [eta2]
    []
    [eta3]
    []

    [grad_eta1]
    []
    [grad_eta2]
    []
    [grad_eta3]
    []
[]

[ICs]
    [eta1_IC]
        type = FunctionIC
        variable = eta1
        function = eta1_IC

    []
    [eta2_IC]
        type = FunctionIC
        variable = eta2
        function = eta2_IC
    []
    [eta3_IC]
        type = FunctionIC
        variable = eta3
        function = eta3_IC
    []
[]

[Functions]
    [line1]
        type = ParsedFunction
        symbol_names = 'theta'
        symbol_values = '${fparse ${pi} / 6 }'
        expression = '(y-${y0})-tan(theta)*(x-${x0})'
    []
    [line2]
        type = ParsedFunction
        symbol_names = 'theta'
        symbol_values = '${fparse 5 * ${pi} / 6 }'
        expression = '(y-${y0})-tan(theta)*(x-${x0})'
    []
    [line3]
        type = ParsedFunction
        symbol_names = 'theta'
        symbol_values = '${fparse -1 * ${pi} / 2 }'
        expression = '(y-${y0})-tan(theta)*(x-${x0})'
    []

    [eta3_IC]
        type = ParsedFunction
        symbol_names = 'l1 l2'
        symbol_values = 'line1 line2'
        expression = '(0.5+0.5*tanh(2*l1/${l}))*(0.5+0.5*tanh(2*l2/${l}))'
    []
    [eta1_IC]
        type = ParsedFunction
        symbol_names = 'l2 l3'
        symbol_values = 'line2 line3'
        expression = '(0.5-0.5*tanh(2*l2/${l}))*(0.5+0.5*tanh(2*l3/${l}))'
    []
    [eta2_IC]
        type = ParsedFunction
        symbol_names = 'l1 l3'
        symbol_values = 'line1 line3'
        expression = '(0.5-0.5*tanh(2*l1/${l}))*(0.5-0.5*tanh(2*l3/${l}))'
    []

[]

[Kernels]
    [dt_eta1]
        type = ADTimeDerivative
        variable = eta1
    []
    [double_well1]
        type = ADAllenCahn
        f_name = F
        variable = eta1
    []
    [gradient_penalty1]
        type = ADACInterface
        variable = eta1
        kappa_name = kappa_int
    []
    [stabilize_moelans_1]
        type = ADACStabilize
        variable = eta1
        v = grad_eta1
        del_kappa_name = 'del_kappa'
        thresh = 1e-4
    []

    [dt_eta2]
        type = ADTimeDerivative
        variable = eta2
    []
    [double_well2]
        type = ADAllenCahn
        f_name = F
        variable = eta2
    []
    [gradient_penalty2]
        type = ADACInterface
        variable = eta2
        kappa_name = kappa_int
    []
    [stabilize_moelans_2]
        type = ADACStabilize
        variable = eta2
        v = grad_eta2
        del_kappa_name = 'del_kappa'
        thresh = 1e-4
    []
    [dt_eta3]
        type = ADTimeDerivative
        variable = eta3
    []
    [double_well3]
        type = ADAllenCahn
        f_name = F
        variable = eta3
    []
    [gradient_penalty3]
        type = ADACInterface
        variable = eta3
        kappa_name = kappa_int
    []
    [stabilize_moelans_3]
        type = ADACStabilize
        variable = eta3
        v = grad_eta3
        del_kappa_name = 'del_kappa'
        thresh = 1e-4
    []

    [grad_eta1_magnitude]
        type = ADGradientMagnitude
        v = eta1
        variable = grad_eta1
    []
    [grad_eta2_magnitude]
        type = ADGradientMagnitude
        v = eta2
        variable = grad_eta2
    []
    [grad_eta3_magnitude]
        type = ADGradientMagnitude
        v = eta3
        variable = grad_eta3
    []
[]

[AuxVariables]
    [bnds]
    []
[]

[AuxKernels]
    [bounds_calc]
        type = BndsCalcAux
        variable = bnds
        v = 'eta1 eta2 eta3'
    []
[]

[Materials]
    [constants]
        type = ADGenericConstantMaterial
        prop_names = 'sigma_int l'
        prop_values = '1.0 ${l}'
    []
    [L]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        property_name = 'L'
        expression = '1.0'
    []

    [f0]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        property_name = 'f0'
        # expression = 'eta1^2*(1-eta1^2) + eta2^2*(1-eta2^2) + eta3^2*(1-eta3^2)'
        expression = 'eta1^4/4 - eta1^2/2 + 1.5*eta1^2*(eta2^2 + eta3^2)/2 +
                      eta2^4/4 - eta2^2/2 + 1.5*eta2^2*(eta1^2 + eta3^2)/2 +
                      eta3^4/4 - eta3^2/2 + 1.5*eta3^2*(eta1^2 + eta2^2)/2 + 1/4'
    []
    [kappa_int]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'sigma_int l'
        property_name = 'kappa_int'
        expression = '3*sigma_int*l/4'
    []
    [mu_int]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'sigma_int l'
        property_name = 'mu_int'
        expression = 'sigma_int*6/l'
        outputs = exodus
    []
    [sigma_s]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'sigma_int l f_chem(eta1,eta2,eta3) mu1:=D[f_chem,eta1] mu2:=D[f_chem,eta2] mu3:=D[f_chem,eta3]'
        property_name = 'sigma_s'
        expression = 'max(max(max(sigma_int,abs(mu1*l/${f_stab})),abs(mu2*l/${f_stab})),abs(mu3*l/${f_stab}))'
        outputs = 'exodus'
    []
    [mu_s]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'sigma_s(eta1,eta2,eta3) l'
        property_name = 'mu_s'
        expression = 'sigma_s*6/l'
    []
    [del_kappa]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'sigma_s(eta1,eta2,eta3) sigma_int l'
        property_name = 'del_kappa'
        expression = '-3*(sigma_s-sigma_int)*l/4'
        outputs = 'exodus'
    []
    [h3]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        property_name = 'h3'
        expression = 'eta3^2/(eta1^2+eta2^2+eta3^2)'
    []
    [f_chem]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'h3(eta1,eta2,eta3) l'
        property_name = f_chem
        expression = 'h3*${V}/l'
        outputs = exodus
    []
    [F]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'eta1 eta2 eta3 grad_eta1 grad_eta2 grad_eta3'
        material_property_names = 'mu_s f0(eta1,eta2,eta3) f_chem(eta1,eta2,eta3) '
        property_name = F
        expression = 'mu_s*f0 + f_chem'
        outputs = exodus
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
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = 'hypre    boomeramg'
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    # petsc_options_value = 'lu superlu_dist'
    automatic_scaling = true
    nl_max_its = 10
    # line_search = none
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1e-1
        iteration_window = 2
        optimal_iterations = 7
        growth_factor = 1.25
        cutback_factor = 0.8
    []
    end_time = ${fparse abs( 1.5 * r0 / ${V} ) }
    # num_steps = 0
[]
[Debug]
    show_var_residual_norms = true
[]

[Outputs]
    exodus = true
    csv = true
[]