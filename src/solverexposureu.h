/*
 * solverexposureu.h
 *
 *  Created on: 25-apr-2020
 *      Author: P. Sashittal
 */

#ifndef SOLVEREXPOSUREU_H
#define SOLVEREXPOSUREU_H

#include <boost/numeric/ublas/matrix.hpp>
#include <gurobi_c++.h>

class SolverExposureU
{
public:
  typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
  
  /// Constructor
  /// @param F Frequency matrix
  /// @param B Genotype matrix
  /// @param nrThreads Number of threads
  /// @param env Gurobi environment
  SolverExposureU(const BoostDoubleMatrix& F,
                  const BoostDoubleMatrix& B,
                  int nrThreads,
                  int nrStrains,
                  GRBEnv env);
  
  BoostDoubleMatrix solve();
  
protected:
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  
  void init();
  
  /// Initialize variables
  void initVariables();
  
  /// Initialize constraints
  void initConstraints();
  
  /// Initialize objective function
  void initObjective();
  
protected:
  /// Frequency matrix (mutations by samples)
  const BoostDoubleMatrix& _F;
  /// Genotype matrix (mutations by strains)
  const BoostDoubleMatrix& _B;
  /// Number of threads
  const int _nrThreads;
  /// Gurobi environment
  GRBEnv _env;
  /// Gurobi model
  GRBModel _model;
  /// _varBU[i][p] is the inferred frequency of mutation i in sample p
  VarMatrix _varBU;
  /// _varU[j][p] is the mixture proportion of strain j in sample p
  VarMatrix _varU;
  /// _k is the maximum number of strains allowed
  const int _k;
  /// _varK is the indicator variable for a strain being used
  VarArray _varK;
};

#endif // SOLVEREXPOSUREU_H
