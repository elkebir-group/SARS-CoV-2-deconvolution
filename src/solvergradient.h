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
  /// Configuration parameters
  struct Param : public Solver::Param
  {
    Param()
      : Solver::Param()
      , _maxIter(100)
      , _epsilon(0.01)
      , _lambda(1.)
    {
    }

    /// Maximum number of iterations
    int _maxIter;
    /// Epsilon
    double _epsilon;
    /// Lambda
    double _lambda;
  };
  
  /// Constructor
  /// @param input Input instance
  /// @param param Parameters
  SolverGradient(const InputInstance& input,
                 const Param& param);
    
  /// Solve
  bool solve();
  
  /// Initialize matrix B
  void initB(const DoubleMatrix& B);
  
  /// Initialize matrix U
  void initU(const DoubleMatrix& U);
  
protected:
  typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
  
  /// Generate random U matrix
  void randomizeU();
  
  /// Generate random B matrix
  void randomizeB();

protected:
  /// Parameters
  const Param& _param;
  /// Genotype matrix (mutations by strains)
  BoostDoubleMatrix _boostB;
  /// Mixture matrix (strains by samples)
  BoostDoubleMatrix _boostU;
};

#endif // SOLVERGRADIENT_H
