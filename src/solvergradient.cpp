/*
 * solvergradient.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include "solvergradient.h"
#include "solvergradientb.h"
#include "solvergradientu.h"
#include "solution.h"
#include <fstream>

SolverGradient::SolverGradient(const InputInstance& input,
                               int nrStrains,
                               int nrRestarts,
                               int nrThreads,
                               double epsilon)
  : Solver(input, nrStrains, nrThreads)
  , _nrRestarts(nrRestarts)
  , _epsilon(epsilon)
  , _boostB(input.getNrMutations(), nrStrains)
  , _boostU(nrStrains, input.getNrSamples())
{
}

bool SolverGradient::solve()
{
  using namespace boost::numeric::ublas;
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  
  BoostDoubleMatrix boostF(nrMutations, nrSamples);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      boostF(i, p) = _input.getVaf(i, p);
    }
  }
  
  for (int r = 0; r < _nrRestarts; ++r)
  {
    randomizeB();
    randomizeU();
    
    double lambda = 1.;
    while (true)
    {
      // solve for U
//      SolverGradientU solverU(boostF, _boostB, _nrThreads, _env);
//      _boostU = solverU.solve();
      
//      for (int j = 0; j < _nrStrains; ++j)
//      {
//        for (int p = 0; p < nrSamples; ++p)
//        {
//          _boostU(j, p) = std::max(1e-10, _boostU(j, p));
//        }
//      }
      
      std::cout << "Frob norm: " << norm_frobenius(boostF - prod(_boostB, _boostU)) << std::endl;
      
      // solve for B
      SolverGradientB solverB(boostF, _boostB, _boostU, lambda);
      _boostB = solverB.solve();
      
//      for (int i = 0; i < nrMutations; ++i)
//      {
//        for (int j = 0; j < _nrStrains; ++j)
//        {
//          _boostB(i, j) = std::max(1e-10, _boostB(i, j));
//        }
//      }
      
      std::cout << "Frob norm: " << norm_frobenius(boostF - prod(_boostB, _boostU)) << std::endl;
      
      // check if
      BoostDoubleMatrix BB = element_prod(_boostB, _boostB);
      BoostDoubleMatrix BB_diff = BB - _boostB;
      BoostDoubleMatrix BB_diff_squared = element_prod(BB_diff, BB_diff);
      double max_diff = norm_inf(BB_diff_squared);
      
            
      std::cout << "max_diff => " << max_diff << std::endl;
      
//      if (max_diff < _epsilon)
      {
        break;
      }
      
      lambda *= 1.1;
    }
  }
  
  // set original matrices
  std::ofstream outTmp("doubleB.tsv");
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      if (j > 0)
        outTmp << "\t";
      outTmp << _boostB(i, j);
      _B[i][j] = _boostB(i, j) > 0.5;
    }
    outTmp << std::endl;
  }
  outTmp.close();
  
  for (int j = 0; j < _nrStrains; ++j)
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
      _F[i][p] = boostBU(i, p);
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
  
  std::uniform_int_distribution<> bit(0, 1);

  const int nrMutations = _input.getNrMutations();

  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      _boostB(i,j) = bit(g_rng);
    }
  }
}

void SolverGradient::randomizeU()
{
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
  
  // https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution
  
  // Symmetric Dirichlet with concentration parameter alpha = 1
//  std::gamma_distribution<> gamma(1, 1);
//
//  const int nrSamples = _input.getNrSamples();
//  for (int j = 0; j < _nrStrains; ++j)
//  {
//    double sum = 0.;
//    for (int p = 0; p < nrSamples; ++p)
//    {
//      _boostU(j, p) = gamma(g_rng);
//      sum += _boostU(j, p);
//    }
//
//    for (int p = 0; p < nrSamples; ++p)
//    {
//      _boostU(j, p) /= sum;
//    }
//  }
}
