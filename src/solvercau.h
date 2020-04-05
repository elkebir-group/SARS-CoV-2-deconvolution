/*
 * solvercau.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERCAU_H
#define SOLVERCAU_H

#include <gurobi_c++.h>
#include "inputinstance.h"

class SolverCaU
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param B Genotype matrix
  /// @param env Gurobi environment
  SolverCaU(const InputInstance& input,
            int nrStrains,
            int nrThreads,
            const BoolMatrix& B,
            GRBEnv env);
  
  /// Solves the U problem
  /// @param F Ouput frequency matrix
  /// @param U Output mixture matrix
  /// @param objValue Output object value
  bool solve(DoubleMatrix& F,
             DoubleMatrix& U,
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
  /// Genotype matrix
  const BoolMatrix& _B;
  
  /// Gurobi model
  GRBModel _model;
  /// _varF[i][p] is the inferred frequency of mutation i in sample p
  VarMatrix _varF;
  /// _varU[j][p] is the mixture proportion of strain j in sample p
  VarMatrix _varU;
};

#endif // SOLVERCAU_H
