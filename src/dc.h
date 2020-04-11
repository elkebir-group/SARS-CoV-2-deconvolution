/*
 * dc.h
 *
 *  Created on: 10-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef DC_H
#define DC_H

#include "inputinstance.h"
#include "solver.h"
#include "solution.h"

/// Divide and conquer class
class DC
{
public:
  /// Constructor
  /// @param input Input instance
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  /// @param sampleLimit Sample limit
  DC(const InputInstance& input,
     int nrThreads,
     int timeLimit,
     int sampleLimit);
  
  void solve();
  
private:
  void solve(const InputInstance& input);
  
  void divide(const InputInstance& input,
              IntSet& samples1,
              IntSet& samples2);
  
  Solution mergeSolutions(const InputInstance& input1,
                          const Solution& sol1,
                          const InputInstance& input2,
                          const Solution& sol2) const;
  
protected:
  /// Input instance
  const InputInstance& _input;
  /// Number of threads
  const int _nrThreads;
  /// Time limit (seconds)
  const int _timeLimit;
  /// Sample limit
  const int _sampleLimit;
  /// Filtered down instance
  InputInstance _filteredInput;
  /// Filtered down input by located
  InputInstance::InputInstanceMap _filteredInputByLocation;
};

#endif // DC_H
