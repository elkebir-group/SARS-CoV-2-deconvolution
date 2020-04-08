/*
 * inputinstance.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "inputinstance.h"

void InputInstance::initSampleLocations()
{
  _sampleLocations.clear();
  for (const std::string& sample : _samples)
  {
    StringVector s;
    boost::split(s, sample, boost::is_any_of("|:"));
    _sampleLocations.push_back(s[1]);
  }
}

InputInstance InputInstance::filterSamplesByLocation(const std::string& location,
                                                     IntVector& newToOld,
                                                     IntVector& oldToNew) const
{
  oldToNew = IntVector(_nrSamples, -1);
  
  for (int p = 0; p < _nrSamples; ++p)
  {
    if (_sampleLocations[p] == location)
    {
      newToOld.push_back(p);
    }
  }
  
  InputInstance newInput;
  
  newInput._alt = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  newInput._ref = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  newInput._vaf = DoubleMatrix(_nrMutations, DoubleVector(_nrSamples, 0));
  newInput._mutDetails = _mutDetails;
  newInput._nrSamples = newToOld.size();
  newInput._nrMutations = _nrMutations;
  
  for (int i = 0; i < _nrMutations; ++i)
  {
    for (int pp = 0; pp < newInput._nrSamples; ++pp)
    {
      int p = newToOld[pp];
      newInput._alt[i][pp] = _alt[i][p];
      newInput._ref[i][pp] = _ref[i][p];
      newInput._vaf[i][pp] = _vaf[i][p];
      if (i == 0)
      {
        newInput._samples.push_back(_samples[p]);
      }
    }
  }
  newInput.initSampleLocations();
  return newInput;
}

InputInstance InputInstance::filterSamples(IntVector& newToOld,
                                           IntVector& oldToNew) const
{
  oldToNew = IntVector(_nrSamples, -1);
  // only retain samples with at least one subclonal mutation
  // and at least 90% of mutation loci with adequate number of reads (>= 5)
  for (int p = 0; p < _nrSamples; ++p)
  {
    int count = 0;
    int adequatePos = 0;
    for (int i = 0; i < _nrMutations; ++i)
    {
      if (getMutationStatus(i, p) == MutSubclonal)
      {
        ++count;
      }
      if (getAltCount(i, p) + getRefCount(i, p) >= 5)
      {
        ++adequatePos;
      }
    }
    
    double threshold = _nrMutations / 2.; //_nrMutations / 10. * 8;
    if (count > 0 && adequatePos >= threshold)
    {
      oldToNew[p] = newToOld.size();
      newToOld.push_back(p);
    }
  }
  
  InputInstance newInput;
  
  newInput._alt = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  newInput._ref = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  newInput._vaf = DoubleMatrix(_nrMutations, DoubleVector(_nrSamples, 0));
  newInput._mutDetails = _mutDetails;
  newInput._nrSamples = newToOld.size();
  newInput._nrMutations = _nrMutations;
  
  for (int i = 0; i < _nrMutations; ++i)
  {
    for (int pp = 0; pp < newInput._nrSamples; ++pp)
    {
      int p = newToOld[pp];
      newInput._alt[i][pp] = _alt[i][p];
      newInput._ref[i][pp] = _ref[i][p];
      newInput._vaf[i][pp] = _vaf[i][p];
      if (i == 0)
      {
        newInput._samples.push_back(_samples[p]);
      }
    }
  }
  newInput.initSampleLocations();
  return newInput;
}

InputInstance InputInstance::filterMutations(IntVector& newToOld,
                                             IntVector& oldToNew) const
{
  oldToNew = IntVector(_nrMutations, -1);
  for (int i = 0; i < _nrMutations; ++i)
  {
    int count = 0;
    int subclonalCount = 0;
    for (int p = 0; p < _nrSamples; ++p)
    {
      switch (getMutationStatus(i, p))
      {
        case MutClonal:
          ++count;
          break;
        case MutSubclonal:
          ++count;
          ++subclonalCount;
        default:
          break;
      }
    }
    
    if (count > 1 && subclonalCount > 0)
    {
      oldToNew[i] = newToOld.size();
      newToOld.push_back(i);
    }
  }
  
  InputInstance newInput;
  for (int i : newToOld)
  {
    newInput._alt.push_back(_alt[i]);
    newInput._ref.push_back(_ref[i]);
    newInput._vaf.push_back(_vaf[i]);
    newInput._mutDetails.push_back(_mutDetails[i]);
  }
  newInput._sampleLocations = _sampleLocations;
  newInput._samples = _samples;
  newInput._nrSamples = _nrSamples;
  newInput._nrMutations = newToOld.size();
  
  return newInput;
}

void InputInstance::read(std::istream& inRef, std::istream& inAlt)
{
  _nrSamples = 0;
  _nrMutations = 0;
  _mutDetails.clear();
  _samples.clear();
  _vaf.clear();
  
  // read ref
  
  // header
  std::string line;
  getline(inRef, line);
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t"));
  
  _nrSamples = s.size() - 10;
  if (_nrSamples < 1)
  {
    throw std::runtime_error("Invalid number of samples.");
  }
  for (int i = 10; i < s.size(); ++i)
  {
    _samples.push_back(s[i]);
  }

  while (inRef.good())
  {
    getline(inRef, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != _nrSamples + 10)
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
    
    _ref.push_back(IntVector(_nrSamples, 0));
    for (int p = 0; p < _nrSamples; ++p)
    {
      _ref.back()[p] = boost::lexical_cast<double>(s[10 + p]);
    }
    _mutDetails.push_back(mutDetails);
    
    _nrMutations++;
  }
  
  // alt
  g_lineNumber = 0;
  getline(inAlt, line);
  boost::split(s, line, boost::is_any_of("\t"));
  
  if (_nrSamples != s.size() - 10)
  {
    throw std::runtime_error("Invalid number of samples.");
  }

  while (inAlt.good())
  {
    getline(inAlt, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != _nrSamples + 10)
    {
      throw std::runtime_error(getLineNumber() + "Invalid number of entries.");
    }
       
    _alt.push_back(IntVector(_nrSamples, 0));
    for (int p = 0; p < _nrSamples; ++p)
    {
      _alt.back()[p] = boost::lexical_cast<double>(s[10 + p]);
    }
  }
  
  _vaf = DoubleMatrix(_nrMutations, DoubleVector(_nrSamples, 0.));
  for (int i = 0; i < _nrMutations; ++i)
  {
    for (int p = 0; p < _nrSamples; ++p)
    {
      _vaf[i][p] = (double)_alt[i][p] / (double)(_alt[i][p] + _ref[i][p]);
    }
  }
  
  initSampleLocations();
}

std::istream& operator>>(std::istream& in, InputInstance& input)
{
  input._nrSamples = 0;
  input._nrMutations = 0;
  input._mutDetails.clear();
  input._samples.clear();
  input._vaf.clear();

  // header
  std::string line;
  getline(in, line);
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t"));
  
  input._nrSamples = s.size() - 10;
  if (input._nrSamples < 1)
  {
    throw std::runtime_error("Invalid number of samples.");
  }
  for (int i = 10; i < s.size(); ++i)
  {
    input._samples.push_back(s[i]);
  }

  while (in.good())
  {
    getline(in, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() != input._nrSamples + 10)
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
    
    input._vaf.push_back(DoubleVector(input._nrSamples, std::nan("")));
    for (int p = 0; p < input._nrSamples; ++p)
    {
      if (s[10 + p] != "NULL" && !s[10 + p].empty())
      {
        input._vaf.back()[p] = boost::lexical_cast<double>(s[10 + p]);
      }
    }
    input._mutDetails.push_back(mutDetails);
    
    input._nrMutations++;
  }
  
  input.initSampleLocations();
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const InputInstance& input)
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
  
  for (const std::string& sample : input._samples)
  {
    out << "\t" << sample;
  }
  
  out << std::endl;
  
  const int n = input._nrMutations;
  const int m = input._nrSamples;
  for (int i = 0; i < n; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = input._mutDetails[i];
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

      const double f_ip = input._vaf[i][p];
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
  
  return out;
}

