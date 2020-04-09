/*
 * solvermiqp.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERMIQP_H
#define SOLVERMIQP_H

#include "solver.h"

class SolverMiqp : public Solver
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  SolverMiqp(const InputInstance& input,
             int nrStrains,
             int nrThreads,
             int timeLimit);
  
  /// Solve
  virtual bool solve();
  
protected:
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  /// Gurobi variable 3D matrix
  typedef std::vector<VarMatrix> Var3Matrix;
  
  /// Initialize
  virtual void init();
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective() = 0;
  
protected:
  /// Time limit (seconds)
  const int _timeLimit;
  /// Gurobi model
  GRBModel _model;
  /// _varF[i][p] is the inferred frequency of mutation i in sample p
  VarMatrix _varF;
  /// _varB[i][j] indicates whether mutation i is present in strain j
  VarMatrix _varB;
  /// _varU[j][p] is the mixture proportion of strain j in sample p
  VarMatrix _varU;
  /// _varBU[i][p][j] = _varB[i][j] * _varU[j][p]
  Var3Matrix _varBU;
};

#endif // SOLVERMIQP_H
