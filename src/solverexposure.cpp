/*
 * solverexposure.cpp
 *
 *  Created on: 25-apr-2020
 *      Author: P. Sashittal
 */

#include <cmath>
#include "solverexposure.h"
#include "solverexposureu.h"
#include "solution.h"
#include <fstream>

SolverExposure::SolverExposure(const InputInstance& input,
                               const Param& param,
                               std::string outputPrefix)
: Solver(input, param)
, _param(param)
, _boostB()
, _boostU()
, _outputPrefix(outputPrefix)
{
}

SolverExposure::SolverExposure(const InputInstance& input,
                               const Param& param,
                               const int nrStrains,
                               std::string outputPrefix)
: Solver(input, param, nrStrains)
, _param(param)
, _boostB()
, _boostU()
, _outputPrefix(outputPrefix)
{
}


void SolverExposure::initB(const DoubleMatrix& B)
{
  const int nrMutations = _input.getNrMutations();
  const int nrStrains = B[0].size();
  
  _boostB = BoostDoubleMatrix(_input.getNrMutations(), nrStrains);
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      _boostB(i,j) = B[i][j];
    }
  }
}

bool SolverExposure::solve()
{
  using namespace boost::numeric::ublas;
  
  const int nrMutations = _input.getNrMutations();
  const int nrSamples = _input.getNrSamples();
  const int nrStrains = _boostB.size2();
  
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
  
  double frobNorm = 0;
  
  // solve for U
  SolverExposureU solverU(boostF, _boostB, _param._nrThreads, _param._nrStrains, _env);
  _boostU = solverU.solve();
  
  frobNorm = norm_frobenius(element_prod(boostM, boostF - prod(_boostB, _boostU)));
  std::cout << "Frob norm 1 after updating U: " << frobNorm << std::endl;
  
  
  std::cout << "nrStrains is " << nrStrains << std::endl;
  std::cout << "size of _B is " << _B.size() << " x " << _B[0].size() << std::endl;
  std::cout << "size of _U is " << _U.size() << " x " << _U[0].size() << std::endl;
  std::cout << "size of boostB is " << _boostB.size1() << " x " << _boostB.size2() << std::endl;
  std::cout << "size of boostU is " << _boostU.size1() << " x " << _boostU.size2() << std::endl;
  
  
  _doubleB = DoubleMatrix(nrMutations, DoubleVector(nrStrains));
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      _doubleB[i][j] = _boostB(i, j);
      _B[i][j] = _boostB(i, j) > 0.5;
    }
  }
  
  for (int j = 0; j < nrStrains; ++j)
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
