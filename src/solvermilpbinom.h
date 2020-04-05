/*
 * solvermilpbinom.h
 *
 *  Created on: 5-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERMILPBINOM_H
#define SOLVERMILPBINOM_H

#include "solvermilp.h"

class SolverMilpBinom : public SolverMilp
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  /// @param nrBreakpoints Number of breakpoints
  SolverMilpBinom(const InputInstance& input,
                  int nrStrains,
                  int nrThreads,
                  int timeLimit,
                  int nrBreakpoints);
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
protected:
  /// Number of segments
  const int _nrBreakpoints;
  /// _varLambda[i][p][l]
  Var3Matrix _varLambda;
  /// Coordinates
  DoubleVector _coord;
  /// Epsilon
  constexpr static double _eps = 1e-3;
};

#endif // SOLVERMILPBINOM_H
