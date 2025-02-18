//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADACStabilize.h"
#include <cmath>
registerMooseObject("PhaseFieldApp", ADACStabilize);

InputParameters
ADACStabilize::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Add the interface stabilization term for stabilized Allen-Cahn.");
  params.addRequiredCoupledVar("v", "Coupled variable with the magnitude of the order parameter gradient");
  params.addParam<Real>("thresh",1e-3,"The threshold for when to turn off stabilization");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addRequiredParam<MaterialPropertyName>("del_kappa_name", "The excess gradient penalty needed for stabilization");
  return params;
}

ADACStabilize::ADACStabilize(const InputParameters & parameters)
  : ADKernel(parameters),
    _v(adCoupledValue("v")),
    _grad_v(coupledGradient("v")),
    _thresh(getParam<Real>("thresh")),
    _prop_L(getADMaterialProperty<Real>("mob_name")),
    _del_kappa(getADMaterialProperty<Real>("del_kappa"))
{
}

ADReal
ADACStabilize::computeQpResidual()
{
    if (MetaPhysicL::raw_value(_v[_qp]) <= _thresh)
        return 0.0;
    return _prop_L[_qp]*_del_kappa[_qp]*_grad_u[_qp]*_grad_v[_qp]*_test[_i][_qp]/_v[_qp];
}
