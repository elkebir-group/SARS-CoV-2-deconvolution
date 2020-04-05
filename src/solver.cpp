/*
 * solver.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "solver.h"
#include "solvercaul1.h"
#include "solvercabl1.h"

Solver::Solver(const InputInstance& input,
               int nrStrains,
               int nrThreads)
  : _input(input)
  , _nrStrains(nrStrains)
  , _nrThreads(nrThreads)
  , _objectiveValue(-std::numeric_limits<double>::max())
  , _B(input.getNrMutations(), BoolVector(_nrStrains, false))
  , _U(_nrStrains, DoubleVector(input.getNrSamples(), 0.))
  , _F(input.getNrMutations(), DoubleVector(input.getNrSamples(), 0.))
  , _env()
{
}

void Solver::writeSolF(std::ostream& out) const
{
  out << "pos" << "\t"
      << "ref" << "\t"
      << "alt" << "\t"
      << "gene" << "\t"
      << "N/S" << "\t"
      << "AA" << "\t"
      << "nSRAsubclonal" << "\t"
      << "nSRAclonal" << "\t"
      << "nSRA" << "\t"
      << "nConsensus";
  
  for (const std::string& sample : _input.getSamples())
  {
    out << "\t" << sample;
  }
  
  out << std::endl;
  
  const int n = _input.getNrMutations();
  const int m = _input.getNrSamples();
  for (int i = 0; i < n; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = _input.getMutationDetails(i);
    out << mutDetails_i._pos << "\t"
        << mutDetails_i._refAllele << "\t"
        << mutDetails_i._altAllele << "\t"
        << mutDetails_i._gene << "\t"
        << mutDetails_i._type << "\t"
        << mutDetails_i._aminoAcidSub << "\t"
        << mutDetails_i._nrSubclonalSamples << "\t"
        << mutDetails_i._nrClonalSamples << "\t"
        << mutDetails_i._nrSamples << "\t"
        << mutDetails_i._nrConsensusSamples;
    
    for (int p = 0; p < m; ++p)
    {
      out << "\t";

      const double f_ip = _F[i][p];
      if (std::isnan(f_ip))
      {
        out << "NULL";
      }
      else
      {
        out << f_ip;
      }
    }
    
    out << std::endl;
  }
}

void Solver::writeSolU(std::ostream& out) const
{
  out << "strain";
  for (const std::string& sample : _input.getSamples())
  {
    out << "\t" << sample;
  }
  out << std::endl;
  
  const int nrSamples = _input.getNrSamples();
  for (int j = 0; j < _nrStrains; ++j)
  {
    out << j;
    for (int p = 0; p < nrSamples; ++p)
    {
      out << "\t" << _U[j][p];
    }
    out << std::endl;
  }
}

void Solver::writeSolB(std::ostream& out) const
{
  out << "pos" << "\t"
      << "ref" << "\t"
      << "alt" << "\t"
      << "gene" << "\t"
      << "N/S" << "\t"
      << "AA" << "\t"
      << "nSRAsubclonal" << "\t"
      << "nSRAclonal" << "\t"
      << "nSRA" << "\t"
      << "nConsensus";
  
  for (int j = 0; j < _nrStrains; ++j)
  {
    out << "\tstrain" << j;
  }
  out << std::endl;
  
  const int nrMutations = _input.getNrMutations();
  for (int i = 0; i < nrMutations; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = _input.getMutationDetails(i);
    out << mutDetails_i._pos << "\t"
        << mutDetails_i._refAllele << "\t"
        << mutDetails_i._altAllele << "\t"
        << mutDetails_i._gene << "\t"
        << mutDetails_i._type << "\t"
        << mutDetails_i._aminoAcidSub << "\t"
        << mutDetails_i._nrSubclonalSamples << "\t"
        << mutDetails_i._nrClonalSamples << "\t"
        << mutDetails_i._nrSamples << "\t"
        << mutDetails_i._nrConsensusSamples;
    
    for (int j = 0; j < _nrStrains; ++j)
    {
      out << "\t" << _B[i][j];
    }
    out << std::endl;
  }
}
