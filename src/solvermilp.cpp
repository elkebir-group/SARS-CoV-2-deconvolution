/*
 * solvermilp.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvermilp.h"

SolverMilp::SolverMilp(const InputInstance& input,
                       int nrStrains,
                       int nrThreads,
                       int timeLimit)
  : Solver(input, nrStrains, nrThreads)
  , _timeLimit(timeLimit)
  , _model(_env)
  , _varF()
  , _varB()
  , _varU()
  , _varBU()
{
}

bool SolverMilp::solve()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  init();
//  _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  if (_nrThreads != -1)
  {
    _model.getEnv().set(GRB_IntParam_Threads, _nrThreads);
  }
  if (_timeLimit > 0)
  {
    _model.getEnv().set(GRB_DoubleParam_TimeLimit, _timeLimit);
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
      _F[i][p] = _varF[i][p].get(GRB_DoubleAttr_X);
    }
  }
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      _B[i][j] = (_varB[i][j].get(GRB_DoubleAttr_X) > 0.4);
    }
  }
  
  for (int j = 0; j < _nrStrains; ++j)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      _U[j][p] = _varU[j][p].get(GRB_DoubleAttr_X);
    }
  }
  
  _objectiveValue = _model.getObjective().getValue();
  
  return true;
}

void SolverMilp::init()
{
  initVariables();
  initConstraints();
  initObjective();
}

void SolverMilp::initVariables()
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
  
  _varBU = Var3Matrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varBU[i] = VarMatrix(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      _varBU[i][p] = VarArray(_nrStrains);
      for (int j = 0; j < _nrStrains; ++j)
      {
        snprintf(buf, 1024, "bu:%d:%d:%d", i, p, j);
        _varBU[i][p][j] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
      }
    }
  }
  
  _model.update();
}

void SolverMilp::initConstraints()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  GRBLinExpr sum;
//  GRBQuadExpr sum;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < _nrStrains; ++j)
      {
        _model.addQConstr(_varBU[i][p][j] == _varB[i][j] * _varU[j][p]);
//        sum += _varB[i][j] * _varU[j][p];
//        _model.addConstr(_varBU[i][p][j] <= _varB[i][j]);
//        _model.addConstr(_varBU[i][p][j] <= _varU[j][p]);
//        _model.addConstr(_varBU[i][p][j] >= _varB[i][j] + _varU[j][p] - 1);
        sum += _varBU[i][p][j];
        
      }
//      _model.addQConstr(_varF[i][p] == sum);
      _model.addConstr(_varF[i][p] == sum);
      sum.clear();
    }
  }

//  GRBLinExpr sum2;
//  GRBLinExpr newSum;
//  for (int j = 0; j < _nrStrains; ++j)
//  {
//    sum2.clear();
//    for (int i = 0; i < nrMutations; ++i)
//    {
//      sum2 += _varB[i][j];
//    }
//    if (j > 0)
//    {
//      _model.addConstr(newSum >= sum2);
//    }
//    newSum = sum2;
//  }
}
