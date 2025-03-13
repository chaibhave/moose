[Models]
    [rom]
        type = LibtorchNeuralNetSurrogate
        inputs = 'forces/c_Cr forces/T'
        outputs = 'state/u_1 state/u_2 state/u_3'
        file_path = 'model_fixed.pt'
        x_mean = '5.0959e-01, 1.4234e+03'
        x_std = '2.9128e-01, 4.9342e+02'
        y_mean = '-73955.1406, -76800.4297, -73335.8203'
        y_std = '40081.9141, 39960.3750, 43353.0977'
    []
[]