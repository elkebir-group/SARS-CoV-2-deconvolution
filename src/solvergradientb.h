/*
 * solvergradientb.h
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLVERGRADIENTB_H
#define SOLVERGRADIENTB_H

#include <cmath>
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

			// vector-to-matrix of B
      BoostDoubleMatrix B(nrMutations, nrStrains, 0);
      for (int ii = 0; ii < _n; ++ii)
      {
        int i = ii / nrStrains;
        int j = ii % nrStrains;
        
        B(i, j) = b[ii];
      }
			
			BoostDoubleMatrix BU = prod(B, _U);
			const int nrSamples = _U.size2();

			// compute the objective function and the gradient
      double fb = 0;
      for (int i = 0; i < nrMutations; ++i)
      {
        for (int j = 0; j < nrStrains; ++j)
        {
					double bij = B(i,j);
          int ii = i * nrStrains + j;
					double sign = ( bij - bij * bij > 0 ) * 2.0 - 1.0;

					if ( std::abs(bij - bij * bij) < 1e-4 ) sign = 0;
					
					// penalty term in the gradient
					grad[ii] = _lambda * sign * ( 1 - 2 * bij );
					
					for (int p = 0; p < nrSamples; ++p)
					{
						if (_F(i,p) != -1)
						{
							// update gradient if F is not nan
							grad[ii] += 2 * BU(i,p) * _U(j,p) - 2 * _F(i,p) * _U(j,p);
							
							// update objective function if F is not nan
							fb += (_F(i,p) - BU(i,p)) * (_F(i,p) - BU(i,p));
						}
					}

					fb += _lambda * std::abs(bij * bij - bij);
        }
      }
			
			/*
			for (int i = 0; i < nrMutations; ++i)
			{
				for (int p = 0; p < nrSamples; ++p)
				{
					if (_F(i,p) != NAN)
					{
						fb += (_F(i,p) - BU(i,p)) * (_F(i,p) - BU(i,p));
					}
				}
			}
			*/
			
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
