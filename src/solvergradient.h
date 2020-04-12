/*
 * solvergradient.h
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERGRADIENT_H
#define SOLVERGRADIENT_H

#include "solver.h"
#include <boost/numeric/ublas/matrix.hpp>

class SolverGradient : public Solver
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrRestarts Number of restarts
  /// @param nrThreads Number of threads
  /// @param epsilon Threshold used for termination
  SolverGradient(const InputInstance& input,
                 int nrStrains,
                 int nrRestarts,
                 int nrThreads,
                 double epsilon);
    
  /// Solve
  bool solve();
  
protected:
  typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
  
  /// Generate random U matrix
  void randomizeU();
  
  /// Generate random B matrix
  void randomizeB();

protected:
  /// Number of restarts
  const int _nrRestarts;
  /// Epsilon
  const double _epsilon;
  /// Genotype matrix (mutations by strains)
  BoostDoubleMatrix _boostB;
  /// Mixture matrix (strains by samples)
  BoostDoubleMatrix _boostU;
};

#endif // SOLVERGRADIENT_H
