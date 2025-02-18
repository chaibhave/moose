//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class ADMultiPhaseBarrierMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ADMultiPhaseBarrierMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// property name
  const MaterialPropertyName _f_name;

  /// function value
  ADMaterialProperty<Real> & _prop_F;

  /// function value derivative
  ADMaterialProperty<Real> & _prop_dFdc;
};
