/*
 * solvergradient.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include <cmath>
#include "solvergradient.h"
#include "solvergradientb.h"
#include "solvergradientu.h"
#include "solution.h"
#include <fstream>

SolverGradient::SolverGradient(const InputInstance& input,
                               const Param& param,
															 std::string outputPrefix)
  : Solver(input, param)
  , _param(param)
  , _boostB()
  , _boostU()
  , _outputPrefix(outputPrefix)
{
}

void SolverGradient::initB(const DoubleMatrix& B)
{
  const int nrMutations = _input.getNrMutations();
  
  _boostB = BoostDoubleMatrix(_input.getNrMutations(), _param._nrStrains);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _param._nrStrains; ++j)
    {
      _boostB(i,j) = B[i][j];
    }
  }
}

void SolverGradient::initU(const DoubleMatrix& U)
{
  const int nrSamples = _input.getNrSamples();
  
  _boostU = BoostDoubleMatrix(_param._nrStrains, _input.getNrSamples());
  for (int j = 0; j < _param._nrStrains; ++j)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      _boostU(j,p) = U[j][p];
    }
  }
}

bool SolverGradient::solve()
{
  using namespace boost::numeric::ublas;
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  BoostDoubleMatrix boostF(nrMutations, nrSamples);
  BoostDoubleMatrix boostM(nrMutations, nrSamples, 1.);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      boostF(i, p) = _input.getVaf(i, p);
      if (boostF(i,p) == -1)
      {
        boostM(i, p) = 0;
      }
    }
  }
  
  BoostDoubleMatrix oldB;
  BoostDoubleMatrix oldU;
  
  double frobNorm = 0;
  
  if (_boostB.size1() == 0 && _boostB.size2() == 0)
  {
    randomizeB();
  }
  if (_boostU.size1() == 0 && _boostU.size2() == 0)
  {
    randomizeU();
  }
  
  double lambda = _param._lambdaInit;
  int idx = 0;
  while (true)
  {
    if (idx > 0)
    {
      _boostB = 0.2 * _boostB + 0.8 * oldB;
      std::cout << "performed mixing" << std::endl;
    }
    
    for (int i = 0; i < nrMutations; ++i)
    {
      for (int j = 0; j < _param._nrStrains; ++j)
      {
        _boostB(i, j) = std::max(1e-4, _boostB(i, j));
      }
    }
    
    // solve for U
    SolverGradientU solverU(boostF, _boostB, _param._nrThreads, _env);
    _boostU = solverU.solve();
    
    for (int j = 0; j < _param._nrStrains; ++j)
    {
      for (int p = 0; p < nrSamples; ++p)
      {
        _boostU(j, p) = std::max(1e-4, _boostU(j, p));
      }
    }
    
    frobNorm = norm_frobenius(element_prod(boostM, boostF - prod(_boostB, _boostU)));
    std::cout << "Frob norm 1 after updating U: " << frobNorm << std::endl;
    
    if (idx > 0)
    {
      _boostU = 0.2 * _boostU + 0.8 * oldU;
    }
    
    // solve for B
    SolverGradientB solverB(boostF, _boostB, _boostU, boostM, lambda, _param._nrThreads);
    _boostB = solverB.solve();
    
    for (int i = 0; i < nrMutations; ++i)
    {
      for (int j = 0; j < _param._nrStrains; ++j)
      {
        _boostB(i, j) = std::max(1e-4, _boostB(i, j));
      }
    }

    for (int j = 0; j < _param._nrStrains; ++j)
    {
      for (int p = 0; p < nrSamples; ++p)
      {
        _boostU(j, p) = std::max(1e-4, _boostU(j, p));
      }
    }
    
    oldB = _boostB;
    oldU = _boostU;
    ++idx;
    
    
    frobNorm = norm_frobenius(element_prod(boostM, boostF - prod(_boostB, _boostU)));
    std::cout << "Frob norm 2 after updating B : " << frobNorm << std::endl;
    
    std::cout << "Iteration number -------- " << idx << "  -----------" << std::endl;
    
    std::cout << "lambda : " << lambda << std::endl;
    
    // check if B is integral
    BoostDoubleMatrix BB = element_prod(_boostB, _boostB);
    BoostDoubleMatrix BB_diff = BB - _boostB;
    BoostDoubleMatrix BB_diff_squared = element_prod(BB_diff, BB_diff);
    double max_diff = norm_inf(BB_diff_squared);
    
    std::cout << "max_diff => " << max_diff << std::endl;

    std::cout << "------------------------normalized = " << frobNorm/norm_frobenius(element_prod(boostM, boostF)) << " ----------------\n";

    if (_param._log && idx % 50 == 0)
    {
      std::ofstream outError(_outputPrefix+"_error_"+std::to_string(idx)+".txt");
      BoostDoubleMatrix E = element_prod(boostM, boostF - prod(_boostB, _boostU));
      for (int i = 0; i < nrMutations; ++i)
      {
        for (int j = 0; j < nrSamples; ++j)
        {
          if (j > 0)
            outError << "\t";
          outError << E(i, j);
        }
        outError << std::endl;
      }
      outError.close();

      std::ofstream outTmp(_outputPrefix+"_doubleB_"+std::to_string(idx)+".txt");
      for (int i = 0; i < nrMutations; ++i)
      {
        for (int j = 0; j < _param._nrStrains; ++j)
        {
          if (j > 0)
            outTmp << "\t";
          outTmp << _boostB(i, j);
        }
        outTmp << std::endl;
      }
      outTmp.close();
    }
    
    if ((max_diff < _param._epsilon && frobNorm < 0.05 * norm_frobenius(element_prod(boostM, boostF))) || idx >= _param._maxIter)
    {
      break;
    }
    
    if (lambda < 4) lambda *= _param._lambda;
    //lambda *= 1.1;
    // set original matrices
  }
  
  _doubleB = DoubleMatrix(nrMutations, DoubleVector(_param._nrStrains));
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _param._nrStrains; ++j)
    {
      _doubleB[i][j] = _boostB(i, j);
      _B[i][j] = _boostB(i, j) > 0.5;
    }
  }
  
  for (int j = 0; j < _param._nrStrains; ++j)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      _U[j][p] = _boostU(j, p);
    }
  }
  
  BoostDoubleMatrix boostBU = prod(_boostB, _boostU);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      _BU[i][p] = boostBU(i, p);
    }
  }
  
  return true;
}

void SolverGradient::randomizeB()
{
  //  Solution sol;
  //
  //  std::ifstream inB("sim_B.tsv");
  //  std::ifstream inU("sim_U.tsv");
  //  sol.readSol(inB, inU);
  //
  //  BoolMatrix B = sol.getB();
  //
  //  const int nrMutations = _input.getNrMutations();
  //
  //  for (int i = 0; i < nrMutations; ++i)
  //  {
  //    for (int j = 0; j < _nrStrains; ++j)
  //    {
  //      _boostB(i,j) = B[i][j];
  //    }
  //  }
  
  _boostB = BoostDoubleMatrix(_input.getNrMutations(), _param._nrStrains);
  std::uniform_real_distribution<> bit(0, 0.01);
  
  const int nrMutations = _input.getNrMutations();
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _param._nrStrains; ++j)
    {
      _boostB(i,j) = bit(g_rng);
    }
  }
}

void SolverGradient::randomizeU()
{
  /*
   Solution sol;
   
   std::ifstream inB("sim_B.tsv");
   std::ifstream inU("sim_U.tsv");
   sol.readSol(inB, inU);
   
   const int nrMutations = _input.getNrMutations();
   const int nrSamples = _input.getNrSamples();
   
   DoubleMatrix U = sol.getU();
   for (int j = 0; j < _nrStrains; ++j)
   {
   for (int p = 0; p < nrSamples; ++p)
   {
   _boostU(j,p) = U[j][p];
   }
   }
   */
  
  _boostU = BoostDoubleMatrix(_param._nrStrains, _input.getNrSamples());
  
  // https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution
  // Symmetric Dirichlet with concentration parameter alpha = 1
  std::gamma_distribution<> gamma(1, 1);
  
  const int nrSamples = _input.getNrSamples();
  for (int j = 0; j < _param._nrStrains; ++j)
  {
    double sum = 0.;
    for (int p = 0; p < nrSamples; ++p)
    {
      _boostU(j, p) = gamma(g_rng);
      sum += _boostU(j, p);
    }
    
    for (int p = 0; p < nrSamples; ++p)
    {
      _boostU(j, p) /= sum;
    }
  }
}
