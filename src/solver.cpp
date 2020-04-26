/*
 * solver.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "solver.h"
#include "solvercaul1.h"
#include "solvercabl1.h"

Solver::Solver(const InputInstance& input,
               const Param& param)
  : _input(input)
  , _param(param)
  , _objectiveValue(-std::numeric_limits<double>::max())
  , _objectiveValueLB(-std::numeric_limits<double>::max())
  , _B(input.getNrMutations(), BoolVector(param._nrStrains, false))
  , _U(param._nrStrains, DoubleVector(input.getNrSamples(), 0.))
  , _BU(input.getNrMutations(), DoubleVector(input.getNrSamples(), 0.))
  , _env()
{
}

Solver::Solver(const InputInstance& input,
							 const Param& param,
							 const int nrStrains)
: _input(input)
, _param(param)
, _objectiveValue(-std::numeric_limits<double>::max())
, _objectiveValueLB(-std::numeric_limits<double>::max())
, _B(input.getNrMutations(), BoolVector(nrStrains, false))
, _U(nrStrains, DoubleVector(input.getNrSamples(), 0.))
, _BU(input.getNrMutations(), DoubleVector(input.getNrSamples(), 0.))
, _env()
{
}
