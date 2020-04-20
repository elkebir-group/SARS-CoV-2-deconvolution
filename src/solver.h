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
  /// Configuration parameters
  struct Param
  {
    Param()
      : _nrStrains(1)
      , _nrThreads(1)
    {
    }

    /// Number of strains
    int _nrStrains;
    /// Number of threads
    int _nrThreads;
  };
  
  /// Constructor
  /// @param input Input instance
  /// @param param Parameters
  Solver(const InputInstance& input,
         const Param& param);
	
	Solver(const InputInstance& input,
				 int nrStrains,
				 int nrThreads)
	: _input(input)
	{
		_param._nrThreads = nrThreads;
		_param._nrStrains = nrStrains;
	}
	
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
  
  /// Return inferred frequency matrix BU
  const DoubleMatrix& getBU() const
  {
    return _BU;
  }
  
protected:
  /// Input
  const InputInstance& _input;
  /// Parameters
	Param _param;
  /// Objective balue
  double _objectiveValue;
  /// Objective value lower bound
  double _objectiveValueLB;
  /// Genotype matrix (mutations by strains)
  BoolMatrix _B;
  /// Mixture matrix (strains by samples)
  DoubleMatrix _U;
  /// Inferred frequences (mutations by samples)
  DoubleMatrix _BU;
  /// Gurobi environment
  GRBEnv _env;
};

#endif // SOLVER_H
