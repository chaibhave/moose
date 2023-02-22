//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DampedNestedSolveMaterial.h"
#include "MatrixTools.h"

registerMooseObject("MooseApp", DampedNestedSolveMaterial);

InputParameters
DampedNestedSolveMaterial::validParams()
{
  InputParameters params = NestedSolveMaterial::validParams();
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
  "max_damping_iters",
  100,
  "Maximum number of times to damp the newton solve per iteration");

  params += NestedSolveMaterial::validParams();
  return params;
}

DampedNestedSolveMaterial::DampedNestedSolveMaterial(const InputParameters & parameters)
  : NestedSolveMaterial(parameters),
    _damping_algorithm(getParam<MooseEnum>("damping_algorithm").getEnum<DampAlgoEnum>()),
    _damping_factor(getParam<Real>("damping_factor")),
    _max_damping_iters(getParam<unsigned int>("max_damping_iters")),
    _condition_name(getParam<MaterialName>("conditions")),
    _old_xi_vals(_num_x)
// _condition_mat_prop(getParam<MaterialName>("conditions"))
{
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
DampedNestedSolveMaterial::initialSetup()
{
  _condition = &getMaterialByName(_condition_name);
  _Ri.resize(_num_x);
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    _Ri[m] = &getMaterialByName(_Ri_names[m]);
  }
}

void
DampedNestedSolveMaterial::update_guess(const NestedSolve::Value<> & guess)
{
  // std::cout << "Old x_i values \n";
  for (unsigned int m = 0; m < guess.size(); ++m)
  {
    _old_xi_vals[m] = (*_prop_xi[m])[_qp];
    // std::cout << _old_xi_vals[m] <<"\t";
  }
  // std::cout << "\n";
  switch (_damping_algorithm)
  {
    // Bounded newton
    case DampAlgoEnum::BOUNDED_DAMP:
    {
      // std::cout << "Damped update\n";

      Real _alpha = 1.0;
      unsigned int i = 0;
      do
      {
        get_damped_guess(guess, _alpha);
        _alpha *= _damping_factor;
        _condition->computePropertiesAtQp(_qp);
        if (++i > _max_damping_iters)
        {
          break;
        }
      } while (!((*_C)[_qp]));
      std::cout << "Number of damping loops = " << i << "\n";
      break;
    }
    case DampAlgoEnum::UNDAMPED:
    {
      Real _alpha = 1.0;
      get_damped_guess(guess, _alpha);
      break;
    }
    default:
      mooseError("Damping_algorithm does not match the given options.");
  }
}

void
DampedNestedSolveMaterial::get_damped_guess(const NestedSolve::Value<> & guess, Real & _alpha)
{
  // std::cout << "Damped update with alpha = " << _alpha << "\n";
  for (unsigned int m = 0; m < _num_x; ++m)
  {
    (*_prop_xi[m])[_qp] = _old_xi_vals[m] + _alpha * (guess(m) - _old_xi_vals[m] );
    // std::cout << (*_prop_xi[m])[_qp] << "\t";
  }
  // std::cout << "\n";
}
