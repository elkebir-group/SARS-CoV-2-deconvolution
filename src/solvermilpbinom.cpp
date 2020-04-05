/*
 * solvermilpbinom.cpp
 *
 *  Created on: 5-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvermilpbinom.h"

SolverMilpBinom::SolverMilpBinom(const InputInstance& input,
                                 int nrStrains,
                                 int nrThreads,
                                 int timeLimit,
                                 int nrBreakpoints)
  : SolverMilp(input, nrStrains, nrThreads, timeLimit)
  , _nrBreakpoints(nrBreakpoints)
  , _varLambda()
  , _coord()
{
}

void SolverMilpBinom::initVariables()
{
  SolverMilp::initVariables();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  char buf[1024];
  
  _varLambda = Var3Matrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varLambda[i] = VarMatrix(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      _varLambda[i][p] = VarArray(_nrBreakpoints);
      for (int l = 0; l < _nrBreakpoints; ++l)
      {
        snprintf(buf, 1024, "lambda:%d:%d:%d", i, p, l);
        _varLambda[i][p][l] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
      }
    }
  }
  
  _model.update();
}

void SolverMilpBinom::initConstraints()
{
  SolverMilp::initConstraints();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  // initialize _coord
  _coord = DoubleVector(_nrBreakpoints, 0);
  const double delta = 1. / (_nrBreakpoints - 1);
  for (int l = 0; l < _nrBreakpoints; ++l)
  {
    _coord[l] = l * delta;
  }
  _coord[0] = _eps;
  _coord[_nrBreakpoints - 1] = 1. - _eps;
  
  GRBLinExpr sum, sum2;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int l = 0; l < _nrBreakpoints; ++l)
      {
        sum2 += _varLambda[i][p][l] * _coord[l];
        sum += _varLambda[i][p][l];
      }
      if (_input.getMutationStatus(i, p) == InputInstance::MutAbsent)
      {
//        _model.addConstr(_varF[i][p] == _eps);
//        _model.addConstr(_varLambda[i][p][0] == 1);
      }
      else if (_input.getMutationStatus(i, p) == InputInstance::MutClonal)
      {
//        _model.addConstr(_varF[i][p] == 1 - _eps);
//        _model.addConstr(_varLambda[i][p][_nrBreakpoints - 1] == 1);
      }
      _model.addConstr(sum == 1);
      _model.addConstr(sum2 == _varF[i][p]);
      sum.clear();
    }
  }
  
  _model.update();
}

void SolverMilpBinom::initObjective()
{
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  GRBLinExpr obj;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      const int alt_ip = _input.getAltCount(i, p);
      const int ref_ip = _input.getRefCount(i, p);
      
      obj += lgamma(alt_ip + ref_ip + 1) - lgamma(alt_ip + 1) - lgamma(ref_ip + 1);
      
      for (int l = 0; l < _nrBreakpoints; ++l)
      {
        obj += _varLambda[i][p][l] * (alt_ip * log(_coord[l]) + ref_ip * log(1-_coord[l]));
      }
    }
  }
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
}
