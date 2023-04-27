//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NestedSolveMaterial.h"
#include "MatrixTools.h"

registerMooseObject("MooseApp", NestedSolveMaterial);

InputParameters
NestedSolveMaterial::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription("Performs a nested newton solve on the residual equations provided. "
                             "Solve can use standard Newton-Raphson or damped newton.");

  params.addRequiredParam<std::vector<MaterialPropertyName>>("xi_names",
                                                             "Material properties to solve for.");
  params.addRequiredParam<std::vector<Real>>("xi_IC",
                                             "Initial values of xi in the same order as xi_names.");

  params.addRequiredParam<std::vector<MaterialName>>(
      "Ri",
      "Residual function material properties for the variables. One function is required per "
      "variable.");
  params += NestedSolve::validParams();
  return params;
}

NestedSolveMaterial::NestedSolveMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _xi_names(getParam<std::vector<MaterialPropertyName>>("xi_names")),
    _num_x(_xi_names.size()),
    _prop_xi(_num_x),
    _xi_old(_num_x),
    _xi_IC(getParam<std::vector<Real>>("xi_IC")),
    _Ri_names(getParam<std::vector<MaterialName>>("Ri")),
    _prop_Ri(_num_x),
    _dRidxi(_num_x),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))
{
  if (!(_xi_names.size() == _Ri_names.size()))
  {
    mooseError("Number of residuals must be equal to number of xi material properties");
  }
  // declare material property variables to solve for
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _prop_xi[m] = &declareProperty<Real>(_xi_names[m]);
    _xi_old[m] = &getMaterialPropertyOld<Real>(_xi_names[m]);
  }

  // Obtain residual expressions
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _prop_Ri[m] = &getMaterialPropertyByName<Real>(_Ri_names[m]);
  }

  // Calculate derivatives of residual functions wrt _prop_xi
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _dRidxi[m].resize(_num_x);
    for (unsigned int n = 0; n < _num_x; ++n)
    {
      _dRidxi[m][n] = &getMaterialPropertyDerivative<Real>(_Ri_names[m], _xi_names[n]);
    }
  }
}

void
NestedSolveMaterial::initQpStatefulProperties()
{
  for (unsigned int m = 0; m < _num_x; ++m)
    (*_prop_xi[m])[_qp] = _xi_IC[m];
}

void
NestedSolveMaterial::initialSetup()
{
  _Ri.resize(_num_x);
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _Ri[m] = &getMaterialByName(_Ri_names[m]);
  }
}

void
NestedSolveMaterial::computeQpProperties()
{
  for (unsigned int m = 0; m < _num_x; ++m)
    (*_prop_xi[m])[_qp] = _xi_IC[m];

  // parameters for nested Newton iteration
  NestedSolve::Value<> solution(_num_x);

  for (unsigned int m = 0; m < _num_x; ++m)
    solution(m) = _xi_IC[m];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian)
  {
    Real _alpha = update_guess(guess);
    for (unsigned int m = 0; m < _num_x; ++m)
    {
      _Ri[m]->computePropertiesAtQp(_qp);
    }

    jacobian.setZero();
    for (unsigned int m = 0; m < _num_x; ++m)
    {
      residual[m] = (*_prop_Ri[m])[_qp];
      for (unsigned int n = 0; n < _num_x; ++n)
      {
        jacobian(m, n) = (*_dRidxi[m][n])[_qp];
      }
    }
    return _alpha;
  };

  _nested_solve.nonlinear(solution, compute);

  // assign solution to ci
  for (unsigned int m = 0; m < _num_x; ++m)
    (*_prop_xi[m])[_qp] = solution[m];

  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _Ri[m]->computePropertiesAtQp(_qp);
  }

  // if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  //   mooseException("Nested Newton iteration did not converge.");
}

Real
NestedSolveMaterial::update_guess(const NestedSolve::Value<> & guess)
{
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    (*_prop_xi[m])[_qp] = guess(m);
  }
  // Output is unscaled
  return 1.0;
}