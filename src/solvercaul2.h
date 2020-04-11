/*
 * solvercaul1.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERCAUL2_H
#define SOLVERCAUL2_H

#include "solvercau.h"
#include <gurobi_c++.h>

class SolverCaUL2 : public SolverCaU
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param B Genotype matrix
  /// @param env Gurobi environment
  SolverCaUL2(const InputInstance& input,
              int nrStrains,
              int nrThreads,
              const BoolMatrix& B,
              GRBEnv env);
  
  /// Destructor
  virtual ~SolverCaUL2()
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
};

#endif // SOLVERCAUL2_H
