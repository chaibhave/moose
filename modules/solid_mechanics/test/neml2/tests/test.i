r = 0.3
l = 0.06
[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 100
    ny = 100
[]

[Variables]
    [c_Cr]
    []
    [T]
        initial_condition = 1000 #K
    []

    [eta_FCC]
    []
    [eta_BCC]
    []
[]

[ICs]
    [c_Cr_IC]
        type = FunctionIC
        function = 'c_IC_func'
        variable = c_Cr
    []
    [FCC]
        type = RandomIC
        variable = eta_FCC
        max = 0.333
        min = 0.0
    []
    [BCC]
        type = RandomIC
        variable = eta_BCC
        max = 0.333
        min = 0.0
    []
[]

[Functions]
    [c_IC_func]
        type = ParsedFunction
        expression = '0.5*(1-tanh(2*(sqrt(x^2+y^2)-${r})/${l}))'
    []
[]

[NEML2]
    input = 'NEML_file.i' # aaaa why does it have hte same extension
    [aaaa]
        model = 'rom'
        moose_input_types = 'VARIABLE VARIABLE'
        moose_inputs = 'c_Cr T'
        neml2_inputs = 'forces/c_Cr forces/T'
        neml2_outputs = 'state/u_1 state/u_2 state/u_3'
        moose_outputs = 'u_1 u_2 u_3'
        moose_output_types = 'MATERIAL MATERIAL MATERIAL'

        neml2_derivatives = 'state/u_1 forces/c_Cr'
        moose_derivatives = 'du_1_dc_Cr'
        moose_derivative_types = 'MATERIAL'

        export_outputs = 'u_1 u_2 u_3 du_1_dc_Cr'
        export_output_targets = 'exodus;exodus;exodus;exodus'
    []
[]

# [Materials]
#   [du_1_dc_Cr]
#     type = ADDerivativeParsedMaterial
#     coupled_variables = 'c_Cr T'
#     material_property_names = 'u_1(c_Cr,T) deriv:=D[u_1,c_Cr]'
#     expression = 'deriv'
#     outputs = 'exodus'
#   []
# []

[Problem]
    kernel_coverage_check = FALSE
[]

[Executioner]
    type = Transient
    num_steps = 1
[]
[Outputs]
    exodus = true
[]