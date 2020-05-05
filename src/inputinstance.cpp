/*
 * inputinstance.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "inputinstance.h"

InputInstance::InputInstance()
  : _nrSamples(0)
  , _nrMutations(0)
  , _mutDetails()
  , _samples()
  , _sampleLocations()
  , _vaf()
  , _ref()
  , _alt()
{
}

InputInstance::InputInstance(DoubleMatrix F,
                             int depth)
  : _nrSamples()
  , _nrMutations()
  , _mutDetails()
  , _samples()
  , _sampleLocations()
  , _vaf(F)
  , _ref()
  , _alt()
{
  _nrMutations = _vaf.size();
  if (_nrMutations > 0)
  {
    _nrSamples = _vaf[0].size();
  }
  else
  {
    _nrSamples = 0;
  }
  
  for (int i = 0; i < _nrMutations; ++i)
  {
    MutationDetails det;
    det._pos = i + 1;
    det._refAllele = 'A';
    det._altAllele = 'T';
    det._gene = "ORF";
    det._type = "S";
    det._aminoAcidSub = "NULL";
    det._nrSubclonalSamples = 0;
    det._nrClonalSamples = 0;
    det._nrSamples = 0;
    det._nrConsensusSamples = 0;
    _mutDetails.push_back(det);
  }
  
  for (int p = 0; p < _nrSamples; ++p)
  {
    _samples.push_back("sample" + std::to_string(p) + "|Earth");
    _sampleLocations.push_back("Earth");
  }
  
  _alt = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  _ref = IntMatrix(_nrMutations, IntVector(_nrSamples, 0));
  for (int i = 0; i < _nrMutations; ++i)
  {
    for (int p = 0; p < _nrSamples; ++p)
    {
      if (_vaf[i][p] != -1)
      {
        _alt[i][p] = depth * _vaf[i][p];
        _ref[i][p] = depth - _alt[i][p];
      }
    }
  }
}

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

BoolMatrix InputInstance::readInitB(std::istream& in)
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
  
  BoolMatrix Bmat;
  
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
    
    Bmat.push_back(BoolVector(nrStrains, false));
    for (int j = 0; j < nrStrains; ++j)
    {
      Bmat.back()[j] = boost::lexical_cast<double>(s[10+j]) >= 0.5;
    }
  }
  
  std::cout << "number of strains in the input B matrix -- " << nrStrains << std::endl;
  
  return Bmat;
}

BoolMatrix InputInstance::blowupBmat()
{
  std::set<BoolVector> strainSet;
  
  for (int p = 0; p < _nrSamples; ++p)
  {
    int subclonalCount = 0;
    IntVector subclonalIndices;
    BoolVector baseStrain;
    
    //std::cout << "clonal mutations:";
    
    for (int i = 0; i < _nrMutations; ++i)
    {
      if (getMutationStatus(i, p) == MutAbsent)
      {
        baseStrain.push_back(false);
      }
      else if (getMutationStatus(i, p) == MutClonal)
      {
        //std::cout << " " << i;
        baseStrain.push_back(true);
      }
      else if (getMutationStatus(i, p) == MutSubclonal)
      {
        ++subclonalCount;
        subclonalIndices.push_back(i);
        baseStrain.push_back(false);
      }
    }
    
    //std::cout << std::endl;
    
    if (subclonalCount > 16)
    {
      throw std::runtime_error("too many mutations in sample " + _samples[p]);
    }
    
    for (int k = 0; k < pow(2.0, subclonalCount); ++k)
    {
      int caseNum = k;
      std::string blowupIndex;
      
      while (caseNum != 0)
      {
        blowupIndex = (caseNum%2 == 0 ? "0" : "1") + blowupIndex;
        caseNum/=2;
      }
      
      blowupIndex.insert(blowupIndex.begin(), subclonalCount - blowupIndex.length(), '0');
      
      /*
       std::cout << blowupIndex << " -- " << subclonalCount << " [";
       
       for (int index : subclonalIndices)
       {
       std::cout << " " << index;
       }
       
       std::cout << " ]\n";
       */
      
      BoolVector currStrain = baseStrain;
      
      for (int idx = 0; idx < subclonalIndices.size(); ++idx)
      {
        if (blowupIndex[idx] == '1')
        {
          currStrain[subclonalIndices[idx]] = true;
          
          //std::cout << subclonalIndices[idx] << " made true" << std::endl;
        }
      }
      
      strainSet.insert(currStrain);
    }
    
    /*
     std::cout << _samples[p] << " -- " << subclonalCount << " -- " << strainSet.size() << " [";
     
     for (int index : subclonalIndices)
     {
     std::cout << " " << index;
     }
     
     std::cout << " ]\n";
     */
  }
  
  BoolMatrix Bmat(_nrMutations);
  
  for (BoolVector strain : strainSet)
  {
    for (int i = 0; i < _nrMutations; ++i)
    {
      Bmat[i].push_back(strain[i]);
    }
  }
  
  std::cout << "number of strains after blowup -- " << strainSet.size() <<std::endl;
  
  return Bmat;
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
      if (_alt[i][p] + _ref[i][p] == 0)
      {
        _vaf[i][p] = -1;
      }
      else
      {
        _vaf[i][p] = (double)_alt[i][p] / (double)(_alt[i][p] + _ref[i][p]);
      }
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

void InputInstance::write(std::ostream& outRef, std::ostream& outAlt) const
{
  writeTemp(_ref, outRef);
  writeTemp(_alt, outAlt);
}

std::ostream& operator<<(std::ostream& out, const InputInstance& input)
{
  input.writeTemp(input._vaf, out);
  return out;
}

InputInstance::InputInstanceMap InputInstance::splitSamplesByLocation() const
{
  InputInstanceMap inputByLocation;
  
  // get locations
  StringSet locations;
  for (const std::string& sample : getSampleLocations())
  {
    locations.insert(sample);
  }
  
  // get location specific input
  for (const std::string& loc : locations)
  {
    IntVector newToOld, oldToNew;
    inputByLocation[loc] = filterSamplesByLocation(loc, newToOld, oldToNew);
  }
  
  return inputByLocation;
}

InputInstance InputInstance::filter() const
{
  InputInstance filteredInput(*this);
  
  // filter
  while (true)
  {
    IntVector newToOldMutations, oldToNewMutations;
    filteredInput = filteredInput.filterMutations(newToOldMutations, oldToNewMutations);
    std::cerr << "Filtered out " << oldToNewMutations.size() - newToOldMutations.size() << " mutation(s) that are present in at most one sample." << std::endl;
    
    int deltaMuts = newToOldMutations.size() - oldToNewMutations.size();
    
    IntVector newToOldSamples, oldToNewSamples;
    filteredInput = filteredInput.filterSamples(newToOldSamples, oldToNewSamples);
    std::cerr << "Filtered out " << oldToNewSamples.size() - newToOldSamples.size() << " sample(s) that do not contain subclonal mutations." << std::endl;
    
    int deltaSamples = newToOldSamples.size() - oldToNewSamples.size();
    if (deltaMuts == 0 && deltaSamples == 0) break;
  }
  
  std::cerr << filteredInput.getNrMutations() << " mutations in " << filteredInput.getNrSamples() << " samples left." << std::endl;
  
  return filteredInput;
}

InputInstance InputInstance::filterSamples(const IntSet& sampleIndices) const
{
  InputInstance newInstance;
  newInstance._nrMutations = _nrMutations;
  newInstance._nrSamples = sampleIndices.size();
  newInstance._mutDetails = _mutDetails;
  newInstance._vaf = DoubleMatrix(newInstance._nrMutations, DoubleVector(newInstance._nrSamples, 0));
  newInstance._ref = IntMatrix(newInstance._nrMutations, IntVector(newInstance._nrSamples, 0));
  newInstance._alt = IntMatrix(newInstance._nrMutations, IntVector(newInstance._nrSamples, 0));
  
  
  int pp = 0;
  for (int p : sampleIndices)
  {
    newInstance._samples.push_back(_samples[p]);
    for (int i = 0; i < newInstance._nrMutations; ++i)
    {
      newInstance._vaf[i][pp] = _vaf[i][p];
      newInstance._ref[i][pp] = _ref[i][p];
      newInstance._alt[i][pp] = _alt[i][p];
    }
    ++pp;
  }
  
  newInstance.initSampleLocations();
  
  return newInstance;
}
