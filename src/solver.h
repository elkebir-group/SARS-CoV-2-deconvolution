/*
 * solver.h
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <gurobi_c++.h>
#include "inputinstance.h"

class Solver
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  Solver(const InputInstance& input,
         int nrStrains,
         int nrThreads);
  
  /// Destructor
  virtual ~Solver()
  {
  }
    
  /// Solve
  virtual bool solve() = 0;
  
  /// Return objective value
  double getObjectiveValue() const
  {
    return _objectiveValue;
  }
  
  /// Return objective value lower bound
  double getObjectiveValueLB() const
  {
    return _objectiveValueLB;
  }
  
  /// Return genotype matrix B
  const BoolMatrix& getB() const
  {
    return _B;
  }
  
  /// Return mixture matrix U
  const DoubleMatrix& getU() const
  {
    return _U;
  }
  
  /// Return frequency matrix F
  const DoubleMatrix& getF() const
  {
    return _F;
  }
  
protected:
  /// Input
  const InputInstance& _input;
  /// Number of strains
	int _nrStrains;
  /// Number of threads
  const int _nrThreads;
  /// Objective balue
  double _objectiveValue;
  /// Objective value lower bound
  double _objectiveValueLB;
  /// Genotype matrix (mutations by strains)
  BoolMatrix _B;
  /// Mixture matrix (strains by samples)
  DoubleMatrix _U;
  /// Inferred frequences (mutations by samples)
  DoubleMatrix _F;
  /// Gurobi environment
  GRBEnv _env;
};

#endif // SOLVER_H
