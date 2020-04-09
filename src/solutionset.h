/*
* solutionset.h
*
*  Created on: 8-apr-2020
*      Author: M. El-Kebir
*/

#ifndef SOLUTIONSET_H
#define SOLUTIONSET_H

#include "inputinstance.h"
#include "solution.h"

class SolutionSet
{
public:
  /// Default constructor
  SolutionSet(const InputInstance& input);
  
  /// Populate solution set
  /// @param prefix Prefix
  /// @param maxNrStrains Maximum number of strains per location
  void populate(const std::string& prefix,
                int maxNrStrains);
  
  /// Return merged solution
  Solution merge() const;
  
  /// Return input instance
  const InputInstance& getInput()
  {
    return _input;
  }
  
private:
  typedef std::map<std::string, InputInstance> InputInstanceMap;
  
  typedef std::map<std::string, std::map<int, Solution>> SolutionMap;
  
private:
  InputInstance _input;
  
  InputInstanceMap _inputByLocation;
  
  SolutionMap _solutionsByLocation;
};

#endif // SOLUTIONSET_H
