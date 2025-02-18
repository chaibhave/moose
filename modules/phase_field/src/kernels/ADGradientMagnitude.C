//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADGradientMagnitude.h"

registerMooseObject("PhaseFieldApp", ADGradientMagnitude);

InputParameters
ADGradientMagnitude::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Set the kernel variable to a specified component of the gradient of a coupled variable.");
  params.addRequiredCoupledVar("v", "Coupled variable to match gradient magnitude of");
  return params;
}

ADGradientMagnitude::ADGradientMagnitude(const InputParameters & parameters)
  : ADKernel(parameters),
    _v_var(coupled("v")),
    _grad_v(coupledGradient("v"))
{
}

ADReal
ADGradientMagnitude::computeQpResidual()
{
   auto _gradient_magnitude = std::sqrt(_grad_v[_qp]*_grad_v[_qp]);
   return (_u[_qp] - _gradient_magnitude) * _test[_i][_qp];
}
