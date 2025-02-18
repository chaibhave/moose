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

  /// Identity of the coupled gradient magnitude variable
  const unsigned int _v_var;

  /// Gradient of the coupled variable
  const VariableGradient & _grad_v;
};
