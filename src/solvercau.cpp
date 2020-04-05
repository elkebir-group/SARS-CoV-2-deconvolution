/*
 * solvercau.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvercau.h"

SolverCaU::SolverCaU(const InputInstance& input,
                     int nrStrains,
                     int nrThreads,
                     const BoolMatrix& B,
                     GRBEnv env)
  : _input(input)
  , _nrStrains(nrStrains)
  , _nrThreads(nrThreads)
  , _B(B)
  , _model(env)
  , _varF()
  , _varU()
{
}

bool SolverCaU::solve(DoubleMatrix& F, DoubleMatrix& U, double& objValue)
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  init();
  _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  if (_nrThreads != -1)
  {
    _model.getEnv().set(GRB_IntParam_Threads, _nrThreads);
  }
  _model.optimize();

  int status = _model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE)
  {
    return false;
  }
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      F[i][p] = _varF[i][p].get(GRB_DoubleAttr_X);
    }
  }
  
  for (int j = 0; j < _nrStrains; ++j)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      U[j][p] = _varU[j][p].get(GRB_DoubleAttr_X);
    }
  }
  
  objValue = _model.getObjective().getValue();
  
  return true;
}

void SolverCaU::init()
{
  initVariables();
  initConstraints();
  initObjective();
}

void SolverCaU::initVariables()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  char buf[1024];
  
  // 1. U problem
  _varF = VarMatrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varF[i] = VarArray(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      snprintf(buf, 1024, "f:%d:%d", i, p);
      _varF[i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _varU = VarMatrix(_nrStrains);
  for (int j = 0; j < _nrStrains; ++j)
  {
    _varU[j] = VarArray(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      snprintf(buf, 1024, "u:%d:%d", j, p);
      _varU[j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _model.update();
}

void SolverCaU::initConstraints()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  GRBLinExpr sum;
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < _nrStrains; ++j)
      {
        sum += _B[i][j] * _varU[j][p];
      }
      _model.addConstr(_varF[i][p] == sum);
      sum.clear();
    }
  }
}
