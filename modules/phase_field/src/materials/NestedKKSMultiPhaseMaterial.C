//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NestedKKSMultiPhaseMaterial.h"

registerMooseObject("PhaseFieldApp", NestedKKSMultiPhaseMaterial);

InputParameters
NestedKKSMultiPhaseMaterial::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription("Solves for the phase-concentrations of the Kim-Kim-Suzuki model "
                             "using the damped nested solve material.");

  params.addRequiredCoupledVar("global_cs", "The interpolated concentrations c, b, etc.");
  params.addRequiredCoupledVar("all_etas", "Order parameters.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Switching functions in the same order as all_etas.");
  params.addRequiredParam<std::vector<MaterialName>>(
      "Fj_material", "Free energy material objects in the same order as all_etas.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. They must have the same order as Fj_names and global_cs, for "
      "example, c1, c2, b1, b2.");
  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");
  params.addCoupledVar("args", "The coupled variables of free energies.");

  // Damped Nested Solve material duplicate code --> FIX
  MooseEnum _damping_algorithm("BOUNDED_DAMP DAMP_FALSE");
  params.addParam<MooseEnum>("damping_algorithm",
                             _damping_algorithm,
                             "Use nested Newton's method with conditional damping (BOUNDED_DAMP) "
                             "or no damping (DAMP_FALSE).");

  params.addParam<MaterialName>("conditions",
                                "C",
                                "Material property that checks bounds and conditions on the "
                                "material properties being solved for.");
  params.addParam<Real>(
      "damping_factor",
      0.8,
      "Number between 0 and 1 that is used to scale the non-linear iteration step size");
  params.addParam<unsigned int>(
      "max_damping_iters", 100, "Maximum number of times to damp the newton solve per iteration");

  params += NestedSolve::validParams();
  return params;
}

NestedKKSMultiPhaseMaterial::NestedKKSMultiPhaseMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _prop_c(coupledValues("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _num_j(coupledComponents("all_etas")),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_j),
    _Fj_names(getParam<std::vector<MaterialName>>("Fj_material")),
    _prop_Fi(_num_j),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_c * _num_j),
    _ci_old(_num_c * _num_j),
    _ci_IC(getParam<std::vector<Real>>("ci_IC")),
    _dFidci(_num_j),
    _d2Fidcidbi(_num_j),
    _args_names(coupledNames("args")),
    _n_args(coupledComponents("args")),
    _dFidarg(_num_j),
    _d2F1dc1darg(_num_c),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _damping_algorithm(getParam<MooseEnum>("damping_algorithm").getEnum<DampAlgoEnum>()),
    _damping_factor(getParam<Real>("damping_factor")),
    _max_damping_iters(getParam<unsigned int>("max_damping_iters")),
    _condition_name(getParam<MaterialName>("conditions")),
    _nested_solve(NestedSolve(parameters)),
    _old_xi_vals(_num_c * _num_j)
{
  // phase concentrations
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
  {
    _ci_old[m] = &getMaterialPropertyOld<Real>(_ci_names[m]);
    _prop_ci[m] = &declareProperty<Real>(_ci_names[m]);
  }

  // free energies
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_Fi[m] = &getMaterialPropertyByName<Real>(_Fj_names[m]);
  }

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    // derivative of free energies wrt phase concentrations
    _dFidci[m].resize(_num_c);
    _d2Fidcidbi[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _dFidci[m][n] = &getMaterialPropertyDerivative<Real>(_Fj_names[m], _ci_names[m + n * _num_j]);

      _d2Fidcidbi[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _d2Fidcidbi[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fj_names[m], _ci_names[m + n * _num_j], _ci_names[m + l * _num_j]);
      }
    }
  }

  // derivative of free energies wrt coupled variables
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _dFidarg[m].resize(_n_args);

    for (unsigned int n = 0; n < _n_args; ++n)
    {
      _dFidarg[m][n] = &getMaterialPropertyDerivative<Real>(_Fj_names[m], _args_names[n]);
    }
  }

  // second derivatives of F1 wrt c1 and other coupled variables
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _d2F1dc1darg[m].resize(_n_args);

    for (unsigned int n = 0; n < _n_args; ++n)
    {
      _d2F1dc1darg[m][n] =
          &getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[m * _num_j], _args_names[n]);
    }
  }

  // Duplicate code --> FIX
  switch (_damping_algorithm)
  {
    // Bounded newton
    case DampAlgoEnum::BOUNDED_DAMP:
    {
      _C = &getMaterialPropertyByName<Real>(_condition_name);
      break;
    }
    case DampAlgoEnum::UNDAMPED:
    {
      std::cout << "Undamped solve running through damped nested solve material\n";
      break;
    }
    default:
      mooseError("Damping_algorithm does not match the given options.");
  }
}

void
NestedKKSMultiPhaseMaterial::initQpStatefulProperties()
{
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    (*_prop_ci[m])[_qp] = _ci_IC[m];
}

void
NestedKKSMultiPhaseMaterial::initialSetup()
{
  _Fj_mat.resize(_num_j);

  for (unsigned int m = 0; m < _num_j; ++m)
    _Fj_mat[m] = &getMaterialByName(_Fj_names[m]);

  _condition = &getMaterialByName(_condition_name);
}

void
NestedKKSMultiPhaseMaterial::computeQpProperties()
{
  // parameters for nested Newton iteration
  NestedSolve::Value<> solution(_num_c * _num_j);

  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    solution(m) = (*_ci_old[m])[_qp];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian)
  {
    Real _alpha = update_guess(guess);

    for (unsigned int m = 0; m < _num_j; ++m)
      _Fj_mat[m]->computePropertiesAtQp(_qp);

    // assign residual functions
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      // if (((*_prop_c[m])[_qp] <= 0.0))
      //   mooseError("Negative values of global concentrations");
  
      for (unsigned int n = 0; n < _num_j - 1; ++n)
        residual(m * _num_j + n) = (*_dFidci[n][m])[_qp] - (*_dFidci[n + 1][m])[_qp];

      residual((m + 1) * _num_j - 1) = -(*_prop_c[m])[_qp];

      for (unsigned int l = 0; l < _num_j; ++l)
        residual((m + 1) * _num_j - 1) += (*_prop_hj[l])[_qp] * (*_prop_ci[m * _num_j + l])[_qp];
    }

    // fill in the non-zero terms in jacobian
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      // equal chemical potential derivative equations
      for (unsigned int n = 0; n < (_num_j - 1); ++n)
      {
        for (unsigned int l = 0; l < _num_c; ++l)
        {
          jacobian(m * _num_j + n, n + l * _num_j) = (*_d2Fidcidbi[n][m][l])[_qp];
          jacobian(m * _num_j + n, n + l * _num_j + 1) = -(*_d2Fidcidbi[n + 1][m][l])[_qp];
        }
      }

      // concentration conservation derivative equations
      for (unsigned int n = 0; n < _num_j; ++n)
        jacobian((m + 1) * _num_j - 1, m * _num_j + n) = (*_prop_hj[n])[_qp];
    }
    return _alpha;
  };

  _nested_solve.nonlinear(solution, compute);

  // assign solution to ci
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    (*_prop_ci[m])[_qp] = solution[m];

  // Update final values for free energies
  for (unsigned int m = 0; m < _num_j; ++m)
    _Fj_mat[m]->computePropertiesAtQp(_qp);
  _condition->computePropertiesAtQp(_qp);
  if (!((*_C)[_qp]))
    std::cout << "Condition failed after nested solve. \n";

  // if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  // //   {
  // //     // auto state = _nested_solve.getState();
  // //     // std::cout << "\n" << state;
  //     mooseError("Nested Newton iteration did not converge.");
  // //   }
}

Real
NestedKKSMultiPhaseMaterial::update_guess(const NestedSolve::Value<> & guess)
{
  for (unsigned int m = 0; m < guess.size(); ++m)
  {
    _old_xi_vals[m] = (*_prop_ci[m])[_qp];
  }
  Real _alpha = 1.0;
  switch (_damping_algorithm)
  {
    // Bounded newton
    case DampAlgoEnum::BOUNDED_DAMP:
    {
      _alpha = 1.0;
      unsigned int damp_ctr = 0;
      do
      {
        if (damp_ctr > 0)
          _alpha *= _damping_factor;

        get_damped_guess(guess, _alpha);
        _condition->computePropertiesAtQp(_qp);
        if (++damp_ctr > _max_damping_iters)
        {
          break;
        }

      } while (!((*_C)[_qp]));
      break;
    }
    case DampAlgoEnum::UNDAMPED:
    {
      _alpha = 1.0;
      get_damped_guess(guess, _alpha);
      break;
    }
    default:
      mooseError("Damping_algorithm does not match the given options.");
  }
  return _alpha;
}

void
NestedKKSMultiPhaseMaterial::get_damped_guess(const NestedSolve::Value<> & guess, Real & _alpha)
{
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
  {
    (*_prop_ci[m])[_qp] = _old_xi_vals[m] + _alpha * (guess(m) - _old_xi_vals[m]);
  }
}
