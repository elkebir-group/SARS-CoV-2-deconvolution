/*
 * solvercaul1.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvercaul2.h"

SolverCaUL2::SolverCaUL2(const InputInstance& input,
                         int nrStrains,
                         int nrThreads,
                         const BoolMatrix& B,
                         GRBEnv env)
: SolverCaU(input, nrStrains, nrThreads, B, env)
{
}

void SolverCaUL2::initVariables()
{
  SolverCaU::initVariables();
  
  _model.update();
}

void SolverCaUL2::initConstraints()
{
  SolverCaU::initConstraints();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  _model.update();
}

void SolverCaUL2::initObjective()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  GRBQuadExpr obj;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      double obs_f_ip = _input.getVaf(i, p);
      if (!std::isnan(obs_f_ip))
      {
        obj += (obs_f_ip - _varF[i][p]) * (obs_f_ip - _varF[i][p]);
      }
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}
