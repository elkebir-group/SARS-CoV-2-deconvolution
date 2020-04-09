/*
 * solvermilpl1.cpp
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvermilpl1.h"

SolverMiqpL2::SolverMiqpL2(const InputInstance& input,
                           int nrStrains,
                           int nrThreads,
                           int timeLimit)
  : SolverMiqp(input, nrStrains, nrThreads, timeLimit)
{
}


void SolverMiqpL2::initConstraints()
{
  SolverMiqp::initConstraints();
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      double obs_f_ip = _input.getVaf(i, p);
      if (!std::isnan(obs_f_ip))
      {
        switch (_input.getMutationStatus(i, p))
        {
          case InputInstance::MutAbsent:
            _model.addConstr(_varF[i][p] == 0);
            break;
          case InputInstance::MutClonal:
            _model.addConstr(_varF[i][p] == 1);
            break;
          case InputInstance::MutSubclonal:
            _model.addConstr(_varF[i][p] >= 0.05);
            break;
        }
      }
    }
  }
  
  _model.update();
}


void SolverMiqpL2::initObjective()
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
