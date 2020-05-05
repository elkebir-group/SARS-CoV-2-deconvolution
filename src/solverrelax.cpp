/*
 * solverca.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "solverrelax.h"
#include "solvercaul1.h"
#include "solvercaul2.h"

SolverRelax::SolverRelax(const InputInstance& input,
                         BoolMatrix B,
                         int nrStrains,
                         int nrThreads,
                         bool L1,
                         double threshold)
: Solver(input, nrStrains, nrThreads)
, _L1(L1)
, _threshold(threshold)
, _Bmat(B)
{
}

bool SolverRelax::solve()
{
  double objValue;
  
  std::cout << "threshold is " << _threshold << std::endl;
  
  int idx = 0;
  while (true)
  {
    
    std::cout << "round " << idx << std::endl;
    
    if (_L1)
    {
      _B = _Bmat;
      SolverCaUL1 solveU(_input, _nrStrains, _nrThreads, _B, _env);
      
      solveU.solve(_BU, _U, objValue);
      
      std::cout << "objective value: " << objValue << std::endl;
      
    }
    else
    {
      _B = _Bmat;
      SolverCaUL2 solveU(_input, _nrStrains, _nrThreads, _B, _env);
      
      solveU.solve(_BU, _U, objValue);
      
      std::cout << "objective value: " << objValue << std::endl;
      
    }
    
    if (! shrinkB())
    {
      break;
    }
    
    ++idx;
  }
  
  return true;
  
}

bool SolverRelax::shrinkB()
{
  int nrMutations = _input.getNrMutations();
  int nrSamples = _input.getNrSamples();
  BoolMatrix newB(nrMutations);
  DoubleMatrix newU;
  IntSet removeIndices;
  
  for (int i = 0; i < _nrStrains; ++i)
  {
    double sum = 0;
    
    for (int j = 0; j < nrSamples; ++j)
    {
      sum += _U[i][j];
    }
    
    //std::cout << "sum for strain " << i << " -- " << sum << std::endl;
    
    if (sum <= _threshold)
    {
      removeIndices.insert(i);
    }
  }
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < _nrStrains; ++j)
    {
      if (removeIndices.find(j) == removeIndices.end())
      {
        newB[i].push_back(_Bmat[i][j]);
      }
    }
  }
  
  for (int j = 0; j < _nrStrains; ++j)
  {
    if (removeIndices.find(j) == removeIndices.end())
    {
      newU.push_back(_U[j]);
    }
  }
  
  std::cout << "number of strains removed -> " << removeIndices.size() << std::endl;
  
  std::cout << "number of starins remain -> " << _nrStrains - removeIndices.size() << std::endl;
  
  if (removeIndices.size() > 0)
  {
    _Bmat = newB;
    _B = _Bmat;
    _U = newU;
    _nrStrains = _nrStrains - removeIndices.size();
    return true;
  }
  else
  {
    return false;
  }
}
