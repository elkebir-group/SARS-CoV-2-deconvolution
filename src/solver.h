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
  
  void writeSolF(std::ostream& out) const;
  
  void writeSolU(std::ostream& out) const;
  
  void writeSolB(std::ostream& out) const;
  
protected:
  /// Input
  const InputInstance& _input;
  /// Number of strains
  const int _nrStrains;
  /// Number of threads
  const int _nrThreads;
  /// Objective Value
  double _objectiveValue;
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