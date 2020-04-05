/*
 * solvermilpl1.h
 *
 *  Created on: 4-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERMILPL1_H
#define SOLVERMILPL1_H

#include "solvermilp.h"

class SolverMilpL1 : public SolverMilp
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrStrains Number of strains
  /// @param nrThreads Number of threads
  SolverMilpL1(const InputInstance& input,
               int nrStrains,
               int nrThreads);
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
protected:
  /// _varG[i][p] = | F[i][p] - varF[i][p] |
  VarMatrix _varG;
};

#endif // SOLVERMILPL1_H
