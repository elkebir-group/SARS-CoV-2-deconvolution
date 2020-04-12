/*
* solution.cpp
*
*  Created on: 8-apr-2020
*      Author: M. El-Kebir
*/

#include <cmath>
#include "solution.h"

Solution::Solution()
  : _B()
  , _U()
  , _F()
  , _mutDetails()
{
}

Solution::Solution(const BoolMatrix& B,
                   const DoubleMatrix& U,
                   const InputInstance::MutationDetailsVector& mutDetails)
  : _B(B)
  , _doubleB()
  , _U(U)
  , _F()
  , _mutDetails(mutDetails)
{
  const int nrSamples = getNrSamples();
  const int nrMutations = getNrMutations();
  const int nrStrains = getNrStrains();
  assert(_mutDetails.size() == _B.size());
  
  _doubleB = DoubleMatrix(nrMutations, DoubleVector(nrStrains, 0));
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      _doubleB[i][j] = 0;
    }
  }
  
  _F = DoubleMatrix(nrMutations, DoubleVector(nrSamples, 0.));
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < nrStrains; ++j)
      {
        _F[i][p] += _B[i][j] * _U[j][p];
      }
    }
  }
}

Solution::Solution(const DoubleMatrix& doubleB,
                   const DoubleMatrix& U,
                   const InputInstance::MutationDetailsVector& mutDetails)
  : _B()
  , _doubleB(doubleB)
  , _U(U)
  , _F()
  , _mutDetails(mutDetails)
{
  const int nrSamples = getNrSamples();
  const int nrMutations = getNrMutations();
  const int nrStrains = getNrStrains();
  assert(_mutDetails.size() == _B.size());
  
  _B = BoolMatrix(nrMutations, BoolVector(nrStrains, false));
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int j = 0; j < nrStrains; ++j)
    {
      _B[i][j] = _doubleB[i][j] >= 0.5;
    }
  }
  
  _F = DoubleMatrix(nrMutations, DoubleVector(nrSamples, 0.));
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < nrStrains; ++j)
      {
        _F[i][p] += _B[i][j] * _U[j][p];
      }
    }
  }
}

BoolVector Solution::getGenotype(int j) const
{
  assert(0 <= j && j < getNrStrains());
  
  BoolVector res;
  
  const int nrMutations = getNrMutations();
  for (int i = 0; i < nrMutations; ++i)
  {
    res.push_back(_B[i][j]);
  }
  
  return res;
}

BoolVector Solution::getExtGenotype(const InputInstance& input, int j) const
{
  assert(0 <= j && j < getNrStrains());
  
  std::map<int, int> posToIdx;
  
  int nrMutationsOrg = input.getNrMutations();
  for (int i = 0; i < nrMutationsOrg; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = input.getMutationDetails(i);
    posToIdx[mutDetails_i._pos] = i;
  }
  
  BoolVector res(nrMutationsOrg, false);
  
  const int nrMutations = getNrMutations();
  for (int i = 0; i < nrMutations; ++i)
  {
    assert(posToIdx.count(_mutDetails[i]._pos) == 1);
    if (_B[i][j])
    {
      res[posToIdx.find(_mutDetails[i]._pos)->second] = true;
    }
  }
  
  return res;
}

int Solution::getNrPresentMutations() const
{
  const int nrMutations = getNrMutations();
  const int nrSamples = getNrSamples();
  
  int presentCount = 0;
  for (int i = 0; i < nrMutations; ++i)
  {
    bool absent = true;
    for (int p = 0; p < nrSamples; ++p)
    {
      if (_F[i][p] > 0.01)
      {
        absent = false;
        break;
      }
    }
    
    if (!absent)
      ++presentCount;
  }
  
  return presentCount;
}

double Solution::getBIC(const InputInstance& input) const
{
  // n*ln(likelihood/n) +k*num_of_strains*ln(n)
  // where n is num of samples times num of mutations
  double k = getNrStrains();
  double n = (getNrSamples() * getNrMutations());
  double dist = computeDistanceL2(input);
  dist /= n;
  
//  double bic = n * log(dist/n) + k * getNrMutations() * log(n);
  double bic = n * log(dist/n) + k * log(n);
  
  return bic;
}

double Solution::computeDistanceL2(const InputInstance& input) const
{
  const int nrSamples = getNrSamples();
  const int nrMutations = getNrMutations();
  
  assert(nrSamples == input.getNrSamples());
  assert(nrMutations == input.getNrMutations());
  
  double dist = 0;
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      double f_ip = _F[i][p];
      double obs_f_ip = input.getVaf(i, p);
      //if (!std::isnan(obs_f_ip))
      if (obs_f_ip != -1)
      {
        double diff = f_ip - obs_f_ip;
        diff *= diff;
        dist += diff;
      }
//
//      switch (input.getMutationStatus(i, p))
//      {
//        case InputInstance::MutClonal:
//          dist += (1 - f_ip) * (1 - f_ip);
//          break;
//        case InputInstance::MutAbsent:
//          dist += (f_ip) * (f_ip);
//          break;
//        case InputInstance::MutSubclonal:
//          {
//            double obs_f_ip = input.getVaf(i, p);
//            if (!std::isnan(obs_f_ip))
//            {
//              double diff = f_ip - obs_f_ip;
//              diff *= diff;
//              dist += diff;
//            }
//          }
//          break;
//      }
    }
  }
  
  return dist;
}

void Solution::readSol(std::istream& inB, std::istream& inU)
{
  readSolU(inU);
  readSolB(inB);
  
  const int nrSamples = getNrSamples();
  const int nrMutations = getNrMutations();
  const int nrStrains = getNrStrains();
  
  _F = DoubleMatrix(nrMutations, DoubleVector(nrSamples, 0.));
  
  for (int i = 0; i < nrMutations; ++i)
  {
    for (int p = 0; p < nrSamples; ++p)
    {
      for (int j = 0; j < nrStrains; ++j)
      {
        _F[i][p] += _B[i][j] * _U[j][p];
      }
    }
  }
}

void Solution::readSolU(std::istream& in)
{
  // header
  std::string line;
  getline(in, line);
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t"));
  
  int nrSamples = s.size() - 1;
  if (nrSamples < 1)
  {
    throw std::runtime_error("Invalid number of samples.");
  }
  
  while (in.good())
  {
    getline(in, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != nrSamples + 1)
    {
      throw std::runtime_error(getLineNumber() + "Invalid number of entries.");
    }
    
    _U.push_back(DoubleVector(nrSamples, 0));
    for (int p = 0; p < nrSamples; ++p)
    {
      _U.back()[p] = boost::lexical_cast<double>(s[1 + p]);
    }
  }
}

void Solution::readSolB(std::istream& in)
{
  // header
  std::string line;
  getline(in, line);
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t"));
  
  int nrStrains = s.size() - 10;
  if (nrStrains < 1)
  {
    throw std::runtime_error("Invalid number of strains.");
  }
  
  while (in.good())
  {
    getline(in, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != nrStrains + 10)
    {
      throw std::runtime_error(getLineNumber() + "Invalid number of entries.");
    }
    
    InputInstance::MutationDetails mutDetails;
    mutDetails._pos = boost::lexical_cast<int>(s[0]);
    mutDetails._refAllele = s[1][0];
    mutDetails._altAllele = s[2][0];
    mutDetails._gene = s[3];
    mutDetails._type = s[4];
    mutDetails._aminoAcidSub = s[5];
    mutDetails._nrSubclonalSamples = boost::lexical_cast<int>(s[6]);
    mutDetails._nrClonalSamples = boost::lexical_cast<int>(s[7]);
    mutDetails._nrSamples = boost::lexical_cast<int>(s[8]);
    mutDetails._nrConsensusSamples = boost::lexical_cast<int>(s[9]);
    
    _doubleB.push_back(DoubleVector(nrStrains, 0.));
    _B.push_back(BoolVector(nrStrains, false));
    for (int j = 0; j < nrStrains; ++j)
    {
      _doubleB.back()[j] = boost::lexical_cast<double>(s[10 + j]);
      _B.back()[j] = _doubleB.back()[j] >= 0.5;
    }
    _mutDetails.push_back(mutDetails);
  }
}

void Solution::writeSolF(const InputInstance& input,
                         std::ostream& out) const
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
  
  for (const std::string& sample : input.getSamples())
  {
    out << "\t" << sample;
  }
  
  out << std::endl;
  
  const int n = input.getNrMutations();
  const int m = input.getNrSamples();
  for (int i = 0; i < n; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = input.getMutationDetails(i);
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

void Solution::writeSolU(const InputInstance& input,
                         std::ostream& out) const
{
  const int nrStrains = getNrStrains();
  
  out << "strain";
  for (const std::string& sample : input.getSamples())
  {
    out << "\t" << sample;
  }
  out << std::endl;
  
  const int nrSamples = input.getNrSamples();
  for (int j = 0; j < nrStrains; ++j)
  {
    out << j;
    for (int p = 0; p < nrSamples; ++p)
    {
      out << "\t" << _U[j][p];
    }
    out << std::endl;
  }
}

void Solution::writeSolB(const InputInstance& input,
                         std::ostream& out) const
{
  const int nrStrains = getNrStrains();
  
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
  
  for (int j = 0; j < nrStrains; ++j)
  {
    out << "\tstrain" << j;
  }
  out << std::endl;
  
  const int nrMutations = input.getNrMutations();
  for (int i = 0; i < nrMutations; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = input.getMutationDetails(i);
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
    
    for (int j = 0; j < nrStrains; ++j)
    {
      out << "\t" << _B[i][j];
    }
    out << std::endl;
  }
}

