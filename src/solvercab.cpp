/*
 * solvercab.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvercab.h"

SolverCaB::SolverCaB(const InputInstance& input,
                     int nrStrains,
                     int nrThreads,
                     const DoubleMatrix& U,
                     GRBEnv env)
  : _input(input)
  , _nrStrains(nrStrains)
  , _nrThreads(nrThreads)
  , _U(U)
  , _model(env)
  , _varF()
  , _varB()
{
}

bool SolverCaB::solve(DoubleMatrix& F, BoolMatrix& B, double& objValue)
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
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      B[i][j] = (_varB[i][j].get(GRB_DoubleAttr_X) > 0.4);
    }
  }
  
  objValue = _model.getObjective().getValue();
  
  return true;
}

void SolverCaB::init()
{
  initVariables();
  initConstraints();
  initObjective();
}

void SolverCaB::initVariables()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  char buf[1024];
  
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
  
  _varB = VarMatrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varB[i] = VarArray(_nrStrains);
    for (int j = 0; j < _nrStrains; ++j)
    {
      snprintf(buf, 1024, "b:%d:%d", i, j);
      _varB[i][j] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
    }
  }
  
  _model.update();
}

void SolverCaB::initConstraints()
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
        sum += _varB[i][j] * _U[j][p];
      }
      _model.addConstr(_varF[i][p] == sum);
      sum.clear();
    }
  }

//  for (int j = 0; j < _nrStrains; ++j)
//  {
//    for (int i = 0; i < nrMutations; ++i)
//    {
//      sum += _varB[i][j];
//    }
//    _model.addConstr(sum >= 1);
//    sum.clear();
//  }
}
