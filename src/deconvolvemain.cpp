/*
 * deconvolvemain.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "inputinstance.h"
#include "solverca.h"
#include "solvermilpl1.h"
#include "solvermilpbinom.h"
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("ca,C", "coordinate ascent")
    ("strains,k", po::value<int>(), "number of strains")
    ("restarts,N", po::value<int>()->default_value(50), "number of restarts")
    ("threads,T", po::value<int>()->default_value(1), "number of threads")
    ("time-limit,l", po::value<int>()->default_value(-1), "time limit in seconds (-1 is unlimited)")
    ("breakpoints,B", po::value<int>()->default_value(50), "number of breakpoints")
    ("seed,s", po::value<int>()->default_value(0), "random number generator seed")
    ("input", po::value<StringVector>(), "input files (ref and alt read counts)")
    ("output,o", po::value<std::string>()->default_value("out"), "output prefix");

  po::positional_options_description p;
  p.add("input", -1);

  try
  {
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << desc << std::endl;
      return 0;
    }

    int nrStrains = vm["strains"].as<int>();
    int nrRestarts = vm["restarts"].as<int>();
    int nrThreads = vm["threads"].as<int>();
    int timeLimit = vm["time-limit"].as<int>();
    int nrBreakpoints = vm["breakpoints"].as<int>();
    std::string inputFilenameRef = vm["input"].as<StringVector>()[0];
    std::string inputFilenameAlt = vm["input"].as<StringVector>()[1];
    std::string outputPrefix = vm["output"].as<std::string>();

    std::ifstream inFileRef(inputFilenameRef.c_str());
    if (!inFileRef.good())
    {
      std::cerr << "Error: could not open '" << inputFilenameRef << "' for reading." << std::endl;
      return 1;
    }

    std::ifstream inFileAlt(inputFilenameAlt.c_str());
    if (!inFileAlt.good())
    {
      std::cerr << "Error: could not open '" << inputFilenameAlt << "' for reading." << std::endl;
      return 1;
    }

    InputInstance input;
    try
    {
      input.read(inFileRef, inFileAlt);
    }
    catch (const std::runtime_error& e)
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }
    inFileRef.close();
    inFileAlt.close();

    InputInstance filteredInput = input;
    while (true)
    {
      IntVector newToOldMutations, oldToNewMutations;
      filteredInput = filteredInput.filterMutations(newToOldMutations, oldToNewMutations);
      std::cerr << "Filtered out " << oldToNewMutations.size() - newToOldMutations.size() << " mutation(s) that are present in at most one sample." << std::endl;

      int deltaMuts = newToOldMutations.size() - oldToNewMutations.size();

      IntVector newToOldSamples, oldToNewSamples;
      filteredInput = filteredInput.filterSamples(newToOldSamples, oldToNewSamples);
      std::cerr << "Filtered out " << oldToNewSamples.size() - newToOldSamples.size() << " sample(s) that do not contain subclonal mutations." << std::endl;
      std::cerr << filteredInput.getNrMutations() << " mutations in " << filteredInput.getNrSamples() << " samples left." << std::endl;

      int deltaSamples = newToOldSamples.size() - oldToNewSamples.size();
      if (deltaMuts == 0 && deltaSamples == 0) break;
    }
    
    int seed = vm["seed"].as<int>();
    g_rng = std::mt19937(seed);


    std::ofstream outFilteredInput(outputPrefix + "_filtered.tsv");
    outFilteredInput << filteredInput;
    outFilteredInput.close();
    
    std::ofstream outLog(outputPrefix + "_log.tsv");
    outLog << "k\tLB\tUB" << std::endl;
    
    for (int k = 1; k <= nrStrains; ++k)
    {
      Solver* pSolve = nullptr;
      if (vm.count("ca"))
      {
        pSolve = new SolverCa(filteredInput, nrStrains, nrRestarts, nrThreads, true, true);
      }
      else
      {
        pSolve = new SolverMilpL1(filteredInput, k, nrThreads, timeLimit);
  //      pSolve = new SolverMilpBinom(filteredInput, nrStrains, nrThreads, timeLimit, nrBreakpoints);
      }
      if (pSolve->solve())
      {
        outLog << k << "\t" << pSolve->getObjectiveValueLB() << "\t" << pSolve->getObjectiveValue() << std::endl;
        
        std::ofstream outF(outputPrefix + "_k" + std::to_string(k) + "_F.tsv");
        pSolve->writeSolF(outF);
        outF.close();
        
        std::ofstream outU(outputPrefix + "_k" + std::to_string(k) + "_U.tsv");
        pSolve->writeSolU(outU);
        outU.close();
        
        std::ofstream outB(outputPrefix + "_k" + std::to_string(k) + "_B.tsv");
        pSolve->writeSolB(outB);
        outB.close();
      }
      else
      {
        outLog << k << "\t" << "-" << "\t" << "-" << std::endl;
      }
      outLog.flush();
      delete pSolve;
    }
    
    outLog.close();
  }
  catch (const boost::program_options::error& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }

  return 0;
}
