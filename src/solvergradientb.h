/*
 * solvergradientb.h
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERGRADIENTB_H
#define SOLVERGRADIENTB_H

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>

class SolverGradientB
{
public:
  typedef boost::numeric::ublas::matrix<double> BoostDoubleMatrix;
  
  SolverGradientB(const BoostDoubleMatrix& F,
                  const BoostDoubleMatrix& B,
                  const BoostDoubleMatrix& U,
                  double lambda);
  
  BoostDoubleMatrix solve() const;
  
protected:
  class Rosenbrock
  {
  private:
    const int _n;
    const BoostDoubleMatrix& _F;
    const BoostDoubleMatrix& _U;
    const double _lambda;
    BoostDoubleMatrix _Ut;
    BoostDoubleMatrix _F_Ut;
    BoostDoubleMatrix _U_Ut;
    
  public:
    Rosenbrock(const BoostDoubleMatrix& F,
               const BoostDoubleMatrix& U,
               double lambda)
      : _n(F.size1() * U.size1())
      , _F(F)
      , _U(U)
      , _lambda(lambda)
      , _Ut(boost::numeric::ublas::trans(_U))
      , _F_Ut(boost::numeric::ublas::prod(_F, _Ut))
      , _U_Ut(boost::numeric::ublas::prod(_U, _Ut))
    {
    }
    
    double operator()(const Eigen::VectorXd& b,
                      Eigen::VectorXd& grad)
    {
      using namespace boost::numeric::ublas;
      
      const int nrMutations = _F.size1();
      const int nrStrains = _U.size1();

      BoostDoubleMatrix B(nrMutations, nrStrains, 0);
      for (int ii = 0; ii < _n; ++ii)
      {
        int i = ii / nrStrains;
        int j = ii % nrStrains;
        
        B(i, j) = b[ii];
      }
      
      BoostDoubleMatrix B_U_Ut = prod(B, _U_Ut);
      BoostDoubleMatrix BB = element_prod(B, B);
      BoostDoubleMatrix BBB = element_prod(BB, B);
      
      for (int i = 0; i < nrMutations; ++i)
      {
        for (int j = 0; j < nrStrains; ++j)
        {
          int ii = i * nrStrains + j;
          grad[ii] = 2 * B_U_Ut(i, j) - 2 * _F_Ut(i, j) + _lambda * (2 * BBB(i, j) + B(i,j) - 3 * BB(i, j));
        }
      }
      
      double fb = pow(norm_frobenius(_F - prod(B, _U)), 2.);
      fb += pow(norm_frobenius(BB - B), 2.);
      
      return fb;
    }
  };
  
private:
  /// Frequency matrix (mutations by samples)
  const BoostDoubleMatrix& _F;
  /// Genotype matrix (mutations by strains)
  const BoostDoubleMatrix& _B;
  /// Mixture matrix (strains by samples)
  const BoostDoubleMatrix& _U;
  /// Lambda
  const double _lambda;
};

#endif // SOLVERGRADIENTB_H
