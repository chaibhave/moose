//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeMaterialInterface.h"
#include "NestedSolve.h"

class NestedSolveMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  NestedSolveMaterial(const InputParameters & parameters);
  /// @brief Overrideable function to make the input material properties change according to nested solver's guess.
  virtual void update_guess(const NestedSolve::Value<> &guess);


  protected:
  virtual void initQpStatefulProperties() override;
  virtual void initialSetup() override;
  virtual void computeQpProperties() override;
  

  //Material property input variables
  const std::vector<MaterialPropertyName> _xi_names;

  //Number of material properties to solve for
  const unsigned int _num_x;


  std::vector<MaterialProperty<Real> *> _prop_xi;
  std::vector<const MaterialProperty<Real> *> _xi_old;
  const std::vector<Real> _xi_IC;

  //Residuals
  const std::vector<MaterialName> _Ri_names;
  std::vector<const MaterialProperty<Real> *> _prop_Ri;
  std::vector<std::vector<const MaterialProperty<Real> *>> _dRidxi;  

  // //conditions for damped solve
  const MaterialPropertyName _condition_mat_prop;
  // const MaterialProperty<Real> & _condition;

  /// Absolute and relative tolerance of nested Newton iteration
  const Real _abs_tol;
  const Real _rel_tol;

  /// Instantiation of the NestedSolve class
  NestedSolve _nested_solve;

  //Instantiation of the MaterialBase class
  std::vector<MaterialBase *> _Ri;

};