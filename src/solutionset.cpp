/*
* solutionset.cpp
*
*  Created on: 8-apr-2020
*      Author: M. El-Kebir
*/

#include "solutionset.h"
#include <fstream>

SolutionSet::SolutionSet(const InputInstance& input)
  : _input(input)
  , _inputByLocation()
  , _solutionsByLocation()
{
  // filter
  while (true)
  {
    IntVector newToOldMutations, oldToNewMutations;
    _input = _input.filterMutations(newToOldMutations, oldToNewMutations);
    std::cerr << "Filtered out " << oldToNewMutations.size() - newToOldMutations.size() << " mutation(s) that are present in at most one sample." << std::endl;

    int deltaMuts = newToOldMutations.size() - oldToNewMutations.size();

    IntVector newToOldSamples, oldToNewSamples;
    _input = _input.filterSamples(newToOldSamples, oldToNewSamples);
    std::cerr << "Filtered out " << oldToNewSamples.size() - newToOldSamples.size() << " sample(s) that do not contain subclonal mutations." << std::endl;
    std::cerr << _input.getNrMutations() << " mutations in " << _input.getNrSamples() << " samples left." << std::endl;

    int deltaSamples = newToOldSamples.size() - oldToNewSamples.size();
    if (deltaMuts == 0 && deltaSamples == 0) break;
  }
  
  // get locations
  StringSet locations;
  for (const std::string& sample : _input.getSampleLocations())
  {
    locations.insert(sample);
  }
  
  // get location specific input
  for (const std::string& loc : locations)
  {
    IntVector newToOld, oldToNew;
    _inputByLocation[loc] = _input.filterSamplesByLocation(loc, newToOld, oldToNew);
  }
}

void SolutionSet::populate(const std::string& prefix,
                           int maxNrStrains)
{
  // get locations
  StringSet locations;
  for (const std::string& sample : _input.getSampleLocations())
  {
    locations.insert(sample);
  }
  
  // get location specific solutions
  _solutionsByLocation.clear();
  for (const std::string& loc : locations)
  {
    for (int k = 1; k <= maxNrStrains; ++k)
    {
      std::string filenameB = prefix + "_" + loc + "_k" + std::to_string(k) + "_B.tsv";
      std::string filenameU = prefix + "_" + loc + "_k" + std::to_string(k) + "_U.tsv";
      
      std::ifstream inU(filenameU.c_str());
      std::ifstream inB(filenameB.c_str());
      
      if (inU.good() && inB.good())
      {
        Solution sol;
        sol.readSol(inB, inU);
        _solutionsByLocation[loc][k] = sol;
        
        assert(_inputByLocation.count(loc) == 1);
        const InputInstance& input = _inputByLocation[loc];
        
        std::cout << loc << "\t" << k << "\t" << sol.computeDistanceL2(input) << "\t" << sol.getBIC(input) << "\t" << sol.getNrPresentMutations() << std::endl;
      }
    }
  }
}

Solution SolutionSet::merge() const
{
  std::map<BoolVector, int> genotypeToStrainIdx;
  int nrMergedStrains = 0;
  int nrMergedSamples = 0;
  
  BoolMatrix B(_input.getNrMutations());
  DoubleMatrix U;
  
  for (const auto& kv : _solutionsByLocation)
  {
    const std::string& loc = kv.first;
    
    int min_k = -1;
    double min_bic = std::numeric_limits<double>::max();
    Solution min_sol;
    for (const auto& kv2 : kv.second)
    {
      const int k = kv2.first;
      const Solution& sol = kv2.second;
      
      double bic = sol.getBIC(_inputByLocation.find(loc)->second);
      if (bic < min_bic)
      {
        min_bic = bic;
        min_k = k;
        min_sol = sol;
      }
    }
    
    // merge
    assert(min_k != -1);
    IntVector strainIndices;
    for (int j = 0; j < min_k; ++j)
    {
      BoolVector extGenotype = min_sol.getExtGenotype(_input, j);
      if (genotypeToStrainIdx.count(extGenotype) == 0)
      {
        genotypeToStrainIdx[extGenotype] = nrMergedStrains;
        U.push_back(DoubleVector(nrMergedSamples, 0.));
        for (int i = 0; i < _input.getNrMutations(); ++i)
        {
          B[i].push_back(extGenotype[i]);
        }
        ++nrMergedStrains;
      }
      strainIndices.push_back(genotypeToStrainIdx[extGenotype]);
    }
      
    // add columns to U corresponding to new samples
    for (int j = 0; j < nrMergedStrains; ++j)
    {
      for (int p = 0; p < min_sol.getNrSamples(); ++p)
      {
        U[j].push_back(0.);
      }
    }
    
    for (int j = 0; j < min_k; ++j)
    {
      int jj = strainIndices[j];
      for (int p = 0; p < min_sol.getNrSamples(); ++p)
      {
        int pp = nrMergedSamples + p;
        U[jj][pp] = min_sol.getU()[j][p];
      }
    }
    
    nrMergedSamples += min_sol.getNrSamples();
  }
  
  return Solution(B, U, _input.getMutationDetails());
}
