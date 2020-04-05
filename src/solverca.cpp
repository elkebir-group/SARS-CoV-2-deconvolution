/*
 * solverca.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "solverca.h"
#include "solvercaul1.h"
#include "solvercabl1.h"

SolverCa::SolverCa(const InputInstance& input,
                   int nrStrains,
                   int nrRestarts,
                   int nrThreads,
                   bool L1,
                   bool randomizeU)
  : Solver(input, nrStrains, nrThreads)
  , _nrRestarts(nrRestarts)
  , _L1(L1)
  , _randomizeU(randomizeU)
{
}

void SolverCa::randomizeU()
{
  // https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution
  
  // Symmetric Dirichlet with concentration parameter alpha = 1
  std::gamma_distribution<> gamma(1, 1);
  
  const int nrSamples = _input.getNrSamples();
  for (int j = 0; j < _nrStrains; ++j)
  {
    double sum = 0.;
    for (int p = 0; p < nrSamples; ++p)
    {
      _U[j][p] = gamma(g_rng);
      sum += _U[j][p];
    }
    
    for (int p = 0; p < nrSamples; ++p)
    {
      _U[j][p] /= sum;
    }
  }
}

void SolverCa::randomizeB()
{
  std::uniform_int_distribution<> bit(0,1);
  
  const int nrMutations = _input.getNrMutations();
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      _B[i][j] = bit(g_rng);
    }
  }
}

bool SolverCa::solve()
{
  DoubleMatrix bestF;
  DoubleMatrix bestU;
  BoolMatrix bestB;
  double bestObj = std::numeric_limits<double>::max();
  
  for (int r = 0; r < _nrRestarts; ++r)
  {
    if (_randomizeU)
    {
      randomizeU();
  
      double delta = std::numeric_limits<double>::max();
      while (delta > 1e-3)
      {
        double newObjValue = 0;
        SolverCaBL1 solveB(_input, _nrStrains, _nrThreads, _U, _env);
        solveB.solve(_F, _B, newObjValue);
  
        SolverCaUL1 solveU(_input, _nrStrains, _nrThreads, _B, _env);
        solveU.solve(_F, _U, newObjValue);
  
        delta = _objectiveValue - newObjValue;
        _objectiveValue = newObjValue;
  
        std::cerr << "Restart " << r << " -- obj value " << _objectiveValue  << " (" << bestObj << ")" << std::endl;
      }
    }
    else
    {
      randomizeB();
      
      double delta = std::numeric_limits<double>::max();
      while (delta > 1e-3)
      {
        double newObjValue = 0;
        SolverCaUL1 solveU(_input, _nrStrains, _nrThreads, _B, _env);
        solveU.solve(_F, _U, newObjValue);
        
        SolverCaBL1 solveB(_input, _nrStrains, _nrThreads, _U, _env);
        solveB.solve(_F, _B, newObjValue);
        
        delta = _objectiveValue - newObjValue;
        _objectiveValue = newObjValue;
        
        std::cerr << "Restart " << r << " -- obj value " << _objectiveValue  << " (" << bestObj << ")" << std::endl;
      }
    }
    
    if (_objectiveValue < bestObj)
    {
      bestF = _F;
      bestU = _U;
      bestB = _B;
      bestObj = _objectiveValue;
    }
  }
  
  _F = bestF;
  _U = bestU;
  _B = bestB;
  _objectiveValue = bestObj;
  
  return true;
}
