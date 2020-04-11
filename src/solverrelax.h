/*
 * solverca.h
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERRELAX_H
#define SOLVERRELAX_H

#include <gurobi_c++.h>
#include "solver.h"

class SolverRelax : public Solver
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrRestarts Number of restarts
  /// @param nrThreads Number of threads
  /// @param L1 Use L1 norm
  /// @param randomizeU Intialize by randomizing U (or B if false)
  SolverRelax(const InputInstance& input,
					 BoolMatrix B,
           int nrStrains,
           int nrThreads,
           bool L1,
					 double threshold);
    
  /// Solve
  bool solve();
  
protected:
  /// Generate random U matrix
  void randomizeU();
  
  /// Generate random B matrix
  void randomizeB();

	/// shrink B matrix
	bool shrinkB();
	
protected:
  /// Use L1 norm
  const bool _L1;
	/// threshold for B
	const double _threshold;
	/// input Bmatrix
	BoolMatrix _Bmat;
};

#endif // SOLVERRELAX_H

