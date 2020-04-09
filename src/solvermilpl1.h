/*
 * solvermilpl1.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERMILPL1_H
#define SOLVERMILPL1_H

#include "solvermilp.h"

class SolverMiqpL2 : public SolverMiqp
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  SolverMiqpL2(const InputInstance& input,
               int nrStrains,
               int nrThreads,
               int timeLimit);
  
protected:	
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
};

#endif // SOLVERMILPL1_H
