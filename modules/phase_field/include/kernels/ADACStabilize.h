//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

class ADACStabilize : public ADKernel
{
public:
  static InputParameters validParams();

  ADACStabilize(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();


  /// Coupled variable
  const ADVariableValue & _v;

  /// Gradient of the coupled gradient magnitude variable
  const VariableGradient & _grad_v;

  /// Threshold
  const Real _thresh;

  /// Mobility
  const ADMaterialProperty<Real> & _prop_L;

  /// Material property containing the driving force for stabilization
  const ADMaterialProperty<Real> & _del_kappa;
};
