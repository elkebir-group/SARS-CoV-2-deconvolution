/*
 * solverexposureu.cpp
 *
 *  Created on: 25-apr-2020
 *      Author: P. Sashittal
 */

#include <cmath>
#include "solverexposureu.h"

SolverExposureU::SolverExposureU(const BoostDoubleMatrix& F,
                                 const BoostDoubleMatrix& B,
                                 int nrThreads,
																 int nrStrainsMax,
                                 GRBEnv env)
  : _F(F)
  , _B(B)
  , _nrThreads(nrThreads)
  , _env(env)
  , _model(_env)
  , _varBU()
  , _varU()
	, _k(nrStrainsMax)
{
}

SolverExposureU::BoostDoubleMatrix SolverExposureU::solve()
{
  const int nrSamples = _F.size2();
  const int nrStrains = _B.size2();
  
  init();
  //_model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  if (_nrThreads != -1)
  {
    _model.getEnv().set(GRB_IntParam_Threads, _nrThreads);
  }
  _model.optimize();
  
  int status = _model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE)
  {
    throw std::runtime_error("Infeasible QP");
  }
  if (_model.get(GRB_IntAttr_SolCount) == 0)
  {
    throw std::runtime_error("Infeasible QP");
  }
  
  BoostDoubleMatrix U(nrStrains, nrSamples);
  
  for (int j = 0; j < nrStrains; ++j)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      U(j,p) = _varU[j][p].get(GRB_DoubleAttr_X);
    }
  }
  
  return U;
}

void SolverExposureU::init()
{
  initVariables();
  initConstraints();
  initObjective();
}

void SolverExposureU::initVariables()
{
  const int nrMutations = _F.size1();
  const int nrSamples = _F.size2();
  const int nrStrains = _B.size2();
  
  char buf[1024];
  
  // 1. U problem
  _varBU = VarMatrix(nrMutations);
  for (int i = 0; i < nrMutations; ++i)
  {
    _varBU[i] = VarArray(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      snprintf(buf, 1024, "f:%d:%d", i, p);
      _varBU[i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _varU = VarMatrix(nrStrains);
  for (int j = 0; j < nrStrains; ++j)
  {
    _varU[j] = VarArray(nrSamples);
    for (int p = 0; p < nrSamples; ++p)
    {
      snprintf(buf, 1024, "u:%d:%d", j, p);
      _varU[j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
	
	_varK = VarArray(nrStrains);
	for (int j = 0; j < nrStrains; ++j)
	{
		snprintf(buf, 1024, "k:%d", j);
		_varK[j] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
	}
	
  _model.update();
}

void SolverExposureU::initConstraints()
{
  const int nrMutations = _F.size1();
  const int nrSamples = _F.size2();
  const int nrStrains = _B.size2();
  
  GRBLinExpr sum;
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < nrStrains; ++j)
      {
        sum += _B(i, j) * _varU[j][p];
      }
      _model.addConstr(_varBU[i][p] == sum);
      sum.clear();
    }
  }
  
  for (int p = 0; p < nrSamples; ++p)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      sum += _varU[j][p];
    }
    _model.addConstr(sum == 1);
    sum.clear();
  }
	
	for (int p = 0; p < nrSamples; ++p)
	{
		for (int j = 0; j < nrStrains; ++j)
		{
			_model.addConstr(_varK[j] >= _varU[j][p]);
		}
	}
	
	for (int j = 0; j < nrStrains; ++j)
	{
		sum += _varK[j];
	}
	_model.addConstr(sum <= _k);
	sum.clear();
	
	for (int p = 0; p < nrSamples; ++p)
	{
		for (int j = 0; j < nrStrains; ++j)
		{
			_model.addConstr(_varU[j][p] >= 0.05 * _varK[j]);
		}
	}
	
  _model.update();
}

void SolverExposureU::initObjective()
{
  const int nrMutations = _F.size1();
  const int nrSamples = _F.size2();
  
  GRBQuadExpr obj;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      if (_F(i,p) != -1)
      {
        double obs_f_ip = _F(i, p);
        obj += (obs_f_ip - _varBU[i][p]) * (obs_f_ip - _varBU[i][p]);
      }
    }
  }
  
  try
  {
    _model.setObjective(obj, GRB_MINIMIZE);
  }
  catch (GRBException e)
  {
    std::cout << "GUROBI OBJ Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    exit(1);
  }
  
  _model.update();
}
