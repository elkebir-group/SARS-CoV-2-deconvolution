/*
 * solvergradientb.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvergradientb.h"
#include <LBFGS.h>

SolverGradientB::SolverGradientB(const BoostDoubleMatrix& F,
                                 const BoostDoubleMatrix& B,
                                 const BoostDoubleMatrix& U,
                                 double lambda)
  : _F(F)
  , _B(B)
  , _U(U)
  , _lambda(lambda)
{
}

SolverGradientB::BoostDoubleMatrix SolverGradientB::solve() const
{
//  using namespace boost::numeric::ublas;
//
//  // U^T
//  BoostDoubleMatrix Ut = trans(_U);
//
//  // F U^T
//  BoostDoubleMatrix F_Ut = prod(_F, Ut);
//
//  // U U^T
//  BoostDoubleMatrix U_Ut = prod(_U, Ut);
//
//  // B U U^T
//  BoostDoubleMatrix B_U_Ut = prod(_B, U_Ut);
//
//  BoostDoubleMatrix BB = element_prod(_B, _B);
//
//  BoostDoubleMatrix B1 = 2 * F_Ut  + _lambda * 3 * BB;
//
//  BoostDoubleMatrix B2 = 2 * B_U_Ut + _lambda * ( 2 * element_prod(_B, BB) + _B);
//
//  BoostDoubleMatrix newB = element_prod(_B, element_div(B1, B2));
//
//  std::cout << "grad => " << norm_inf(B2 - B1) << std::endl;
  
  // Set up parameters
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 1000;
  
  // Create solver and function object
  LBFGSpp::LBFGSSolver<double> solver(param);
  Rosenbrock fun(_F, _U, _lambda);
  
  const int nrMutations = _F.size1();
  const int nrStrains = _U.size1();
  Eigen::VectorXd b(nrMutations * nrStrains);
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      int ii = i * nrStrains + j;
      b[ii] = _B(i, j);
    }
  }
  
  double fb;
  int niter = solver.minimize(fun, b, fb);
  std::cout << niter << " iterations" << std::endl;
  std::cout << "f(b) = " << fb << std::endl;
  
  BoostDoubleMatrix newB(nrMutations, nrStrains, 0);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      int ii = i * nrStrains + j;
    
      newB(i, j) = b[ii];
    }
  }
  
  /*
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      double b_ij = _B(i,j);
      
      double denom = B_U_Ut(i,j) + _lambda * ( 2 * pow(b_ij, 3.) + b_ij);
      
      double alpha = b_ij / denom;
      
      newB(i, j) = b_ij - alpha * (denom - F_Ut(i,j) - 3 * _lambda * pow(b_ij, 2.));
    }
  }
  */
  return newB;
}
