/*
 * solvermilpl1.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvermilpl1.h"

SolverMilpL1::SolverMilpL1(const InputInstance& input,
                           int nrStrains,
                           int nrThreads,
                           int timeLimit)
  : SolverMilp(input, nrStrains, nrThreads, timeLimit)
{
}

void SolverMilpL1::initVariables()
{
  SolverMilp::initVariables();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  char buf[1024];
  
//  _varG = VarMatrix(nrMutations);
//  for (int i = 0; i < nrMutations; ++i)
//  {
//    _varG[i] = VarArray(nrSamples);
//    for (int p = 0; p < nrSamples; ++p)
//    {
//      snprintf(buf, 1024, "g:%d:%d", i, p);
//      _varG[i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
//    }
//  }
  
  _model.update();
}

void SolverMilpL1::initConstraints()
{
  SolverMilp::initConstraints();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      double obs_f_ip = _input.getVaf(i, p);
      if (!std::isnan(obs_f_ip))
      {
//        _model.addConstr(_varG[i][p] >= obs_f_ip - _varF[i][p]);
//        _model.addConstr(_varG[i][p] >= _varF[i][p] -  obs_f_ip);
        
        if (_input.getMutationStatus(i, p) == InputInstance::MutAbsent)
        {
          _model.addConstr(_varF[i][p] == 0);
//          _model.addConstr(_varG[i][p] == obs_f_ip);
        }
        else if (_input.getMutationStatus(i, p) == InputInstance::MutClonal)
        {
          _model.addConstr(_varF[i][p] == 1);
//          _model.addConstr(_varG[i][p] == 1 - obs_f_ip);
        }
//				else if (_input.getMutationStatus(i, p) == InputInstance::MutSubclonal)
//				{
//					_model.addConstr(_varF[i][p] >= 0.05);
//					//_model.addConstr(_varF[i][p] <= 0.95);
//				}
      }
    }
  }
	
//	GRBLinExpr sum;
//
// for (int i = 0; i < nrMutations; ++i)
// {
//		for (int p = 0; p < nrSamples; ++p)
//		{
//			sum += _varF[i][p];
//		}
//	 _model.addConstr(sum >= 0.05);
//	 sum.clear();
//	}
  
  _model.update();
}


void SolverMilpL1::initObjective()
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
