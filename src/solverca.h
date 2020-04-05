/*
 * solverca.h
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERCA_H
#define SOLVERCA_H

#include <gurobi_c++.h>
#include "solver.h"

class SolverCa : public Solver
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrRestarts Number of restarts
  /// @param nrThreads Number of threads
  /// @param L1 Use L1 norm
  /// @param randomizeU Intialize by randomizing U (or B if false)
  SolverCa(const InputInstance& input,
           int nrStrains,
           int nrRestarts,
           int nrThreads,
           bool L1,
           bool randomizeU);
    
  /// Solve
  bool solve();
  
protected:
  /// Generate random U matrix
  void randomizeU();
  
  /// Generate random B matrix
  void randomizeB();
  
protected:
  /// Number of restarts
  const int _nrRestarts;
  /// Use L1 norm
  const bool _L1;
  /// Randomize U
  const bool _randomizeU;
};

#endif // SOLVERCA_H

