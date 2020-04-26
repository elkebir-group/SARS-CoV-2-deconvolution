/*
 * solverexposure.h
 *
 *  Created on: 25-apr-2020
 *      Author: P. Sashittal
 */

#ifndef SOLVEREXPOSURE_H
#define SOLVEREXPOSURE_H

#include "solver.h"
#include <boost/numeric/ublas/matrix.hpp>

class SolverExposure : public Solver
{
public:
  /// Configuration parameters
  struct Param : public Solver::Param
  {
    Param()
      : Solver::Param()
    {
    }
  };
  
  /// Constructor
  /// @param input Input instance
  /// @param param Parameters
  SolverExposure(const InputInstance& input,
                 const Param& param,
								 std::string outputPrefix);

	SolverExposure(const InputInstance& input,
								 const Param& param,
								 const int nrStrains,
								 std::string outputPrefix);
	
  /// Solve
  bool solve();
  
  /// Initialize matrix B
  void initB(const DoubleMatrix& B);
	
  const DoubleMatrix& getDoubleB() const
  {
    return _doubleB;
  }
  
protected:
  typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;

protected:
  /// Parameters
  const Param& _param;
  /// Genotype matrix (mutations by strains)
  BoostDoubleMatrix _boostB;
  /// Mixture matrix (strains by samples)
  BoostDoubleMatrix _boostU;
  /// Genotype matrix (double)
  DoubleMatrix _doubleB;
	/// output prefix
	std::string _outputPrefix;
	
};

#endif // SOLVERGRADIENT_H
