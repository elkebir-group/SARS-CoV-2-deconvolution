/*
 * solvercab.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERCAB_H
#define SOLVERCAB_H

#include <gurobi_c++.h>
#include "inputinstance.h"

/// Coordinate ascent solver -- B problem
class SolverCaB
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param U Mixture matrix
  /// @param env Gurobi environment
  SolverCaB(const InputInstance& input,
            int nrStrains,
            int nrThreads,
            const DoubleMatrix& U,
            GRBEnv env);
  
  /// Destructor
  virtual ~SolverCaB()
  {
  };
  
  /// Solves the U problem
  /// @param F Ouput frequency matrix
  /// @param B Output genotype matrix
  /// @param objValue Output objective value
  bool solve(DoubleMatrix& F,
             BoolMatrix& B,
             double& objValue);
  
protected:
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  
  /// Initialize
  virtual void init();
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective() = 0;
  
protected:
  /// Input instance
  const InputInstance& _input;
  /// Number of strains
  const int _nrStrains;
  /// Number of threads
  const int _nrThreads;
  /// Mixture matrix
  const DoubleMatrix& _U;
  
  /// Gurobi model
  GRBModel _model;
  /// _varF[i][p] is the inferred frequency of mutation i in sample p
  VarMatrix _varF;
  /// _varB[i][j] indicates whether mutation i is present in strain j
  VarMatrix _varB;
};

#endif // SOLVERCAB_H
