/*
 * solvercabl1.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERCABL1_H
#define SOLVERCABL1_H

#include "solvercab.h"
#include <gurobi_c++.h>

class SolverCaBL1 : public SolverCaB
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param U Mixture matrix
  /// @param env Gurobi environment
  SolverCaBL1(const InputInstance& input,
              int nrStrains,
              int nrThreads,
              const DoubleMatrix& U,
              GRBEnv env);
  
  /// Destructor
  virtual ~SolverCaBL1()
  {
  };
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
protected:
  /// _varG[i][p] = | F[i][p] - varF[i][p] |
  VarMatrix _varG;
};

#endif // SOLVERCABL1_H
