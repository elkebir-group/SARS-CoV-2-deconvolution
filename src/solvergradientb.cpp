/*
 * solvergradientb.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvergradientb.h"
#include <LBFGS.h>
//#include <LBFGSB.h>

SolverGradientB::SolverGradientB(const BoostDoubleMatrix& F,
                                 const BoostDoubleMatrix& B,
                                 const BoostDoubleMatrix& U,
                                 const BoostDoubleMatrix& M,
                                 double lambda)
  : _F(F)
  , _B(B)
  , _U(U)
  , _M(M)
  , _lambda(lambda)
{
}

SolverGradientB::BoostDoubleMatrix SolverGradientB::solve() const
{
  // Set up parameters
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-8;
  param.max_iterations = 200;
  param.m = 5;
  
 
  // Create solver and function object
  LBFGSpp::LBFGSSolver<double> solver(param);
  Rosenbrock fun(_F, _U, _M, _lambda);
  
  const int nrMutations = _F.size1();
  const int nrStrains = _U.size1();
  Eigen::VectorXd b(nrMutations * nrStrains);
  
  //Eigen::VectorXd ub(nrMutations * nrStrains);
  //Eigen::VectorXd lb(nrMutations * nrStrains); 

  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      int ii = i * nrStrains + j;
      b[ii] = _B(i, j);
      //ub[ii] = 2;
      //lb[ii] = -1;
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
