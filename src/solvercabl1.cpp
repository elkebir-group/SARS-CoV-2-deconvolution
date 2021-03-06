/*
 * solvercabl1.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvercabl1.h"

SolverCaBL1::SolverCaBL1(const InputInstance& input,
                         int nrStrains,
                         int nrThreads,
                         const DoubleMatrix& U,
                         GRBEnv env)
  : SolverCaB(input, nrStrains, nrThreads, U, env)
{
}

void SolverCaBL1::initVariables()
{
  SolverCaB::initVariables();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  char buf[1024];
  
  _varG = VarMatrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varG[i] = VarArray(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      snprintf(buf, 1024, "g:%d:%d", i, p);
      _varG[i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _model.update();
}

void SolverCaBL1::initConstraints()
{
  SolverCaB::initConstraints();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      double obs_f_ip = _input.getVaf(i, p);
      if (!std::isnan(obs_f_ip))
      {
        _model.addConstr(_varG[i][p] >= obs_f_ip - _varF[i][p]);
        _model.addConstr(_varG[i][p] >= _varF[i][p] -  obs_f_ip);
      }
    }
  }
  
  _model.update();
}

void SolverCaBL1::initObjective()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  GRBLinExpr obj;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      obj += _varG[i][p];
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}
