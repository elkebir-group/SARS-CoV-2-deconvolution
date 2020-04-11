/*
 * dc.cpp
 *
 *  Created on: 10-apr-2020
 *      Author: M. El-Kebir
 */

#include "dkm/dkm.hpp"
#include "dc.h"
#include "solvermilpl1.h"
#include "solution.h"
#include "solutionset.h"

DC::DC(const InputInstance& input,
       int nrThreads,
       int timeLimit,
       int sampleLimit)
  : _input(input)
  , _nrThreads(nrThreads)
  , _timeLimit(timeLimit)
  , _sampleLimit(sampleLimit)
  , _filteredInput(_input.filter())
  , _filteredInputByLocation(_filteredInput.splitSamplesByLocation())
{
}

void DC::solve()
{
  for (const auto& kv : _filteredInputByLocation)
  {
    const std::string& loc = kv.first;
    const InputInstance& input = kv.second;
    
    solve(input);
  }
}

void DC::divide(const InputInstance& input,
                IntSet& samples1,
                IntSet& samples2)
{
  const int MUTATION_LIMIT_UB = 30000;
  
  std::uniform_int_distribution<> unif(0, 100);
  
  const int nrSamples = input.getNrSamples();
  const int nrMutations = input.getNrMutations();
  
  assert(nrMutations < MUTATION_LIMIT_UB);
  
  // divide by 2-means
  std::vector<std::array<double, MUTATION_LIMIT_UB>> data;
  for (int p = 0; p < nrSamples; ++p)
  {
    data.push_back(std::array<double, MUTATION_LIMIT_UB>());
    for (int i = 0; i < nrMutations; ++i)
    {
      double val;
      switch (input.getMutationStatus(i, p))
      {
        case InputInstance::MutClonal:
          val = 1.;
          break;
        case InputInstance::MutSubclonal:
          val = input.getVaf(i, p);
          break;
        case InputInstance::MutAbsent:
          val = 0.;
          break;
      }
      data.back()[i] = val;
    }
  }
  
  auto z = std::get<1>(dkm::kmeans_lloyd(data, 2, unif(g_rng)));
  for (int p = 0; p < nrSamples; ++p)
  {
    assert(0 <= z[p] && z[p] <= 1);
    if (z[p] == 0)
    {
      samples1.insert(p);
    }
    else
    {
      samples2.insert(p);
    }
  }
}

Solution DC::mergeSolutions(const InputInstance& input1,
                            const Solution& sol1,
                            const InputInstance& input2,
                            const Solution& sol2) const
{
  Solution mergedSolution;
  return mergedSolution;
}

void DC::solve(const InputInstance& input)
{
  if (input.getNrSamples() > _sampleLimit)
  {
    // divide by 2-means
    IntSet samples1, samples2;
    divide(input, samples1, samples2);
    
    InputInstance input1;
    input1.filterSamples(samples1);
    
    InputInstance input2;
    input2.filterSamples(samples2);
    
    solve(input1);
    solve(input2);
    
//    Solution mergedSol = mergeSolutions(input1, sol1, input2, sol2);
  }
  else
  {
    // keep solving until BIC increases
    int k = 1;
    double min_bic = std::numeric_limits<double>::max();
    Solution best_sol;
    while (true)
    {
      SolverMiqpL2 solver(input, k, _nrThreads, _timeLimit);
      if (solver.solve())
      {
        Solution sol(solver.getB(), solver.getU(), input.getMutationDetails());
        double bic = sol.getBIC(input);
        if (bic < min_bic)
        {
          min_bic = bic;
          best_sol = sol;
          ++k;
        }
        else
        {
          --k;
          break;
        }
      }
    }
  }
}
