//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "NestedSolveMaterial.h"
#include "MatrixTools.h"

class DampedNestedSolveMaterial : public NestedSolveMaterial
{
public:
  static InputParameters validParams();

  DampedNestedSolveMaterial(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  virtual Real update_guess(const NestedSolve::Value<> & guess) override;
  virtual void get_damped_guess(const NestedSolve::Value<> & guess, Real & _alpha);

  // MooseEnum _damping_algorithm;
  enum class DampAlgoEnum
  {
    BOUNDED_DAMP,
    UNDAMPED
  } _damping_algorithm;

  Real _damping_factor;
  const unsigned int _max_damping_iters;

  MaterialName _condition_name;
  MaterialBase * _condition;
  const MaterialProperty<Real> * _C;
  std::vector<Real> _old_xi_vals;
};