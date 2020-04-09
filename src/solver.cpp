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
               int nrStrains,
               int nrThreads)
  : _input(input)
  , _nrStrains(nrStrains)
  , _nrThreads(nrThreads)
  , _objectiveValue(-std::numeric_limits<double>::max())
  , _objectiveValueLB(-std::numeric_limits<double>::max())
  , _B(input.getNrMutations(), BoolVector(_nrStrains, false))
  , _U(_nrStrains, DoubleVector(input.getNrSamples(), 0.))
  , _F(input.getNrMutations(), DoubleVector(input.getNrSamples(), 0.))
  , _env()
{
}
