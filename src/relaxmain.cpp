/*
 * relaxmain.cpp
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "inputinstance.h"
#include "solution.h"
//#include "solverlpl1.h"
#include "solverrelax.h"
#include "solvercaul1.h"
#include "solvercaul2.h"
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("threads,T", po::value<int>()->default_value(1), "number of threads")
    ("input", po::value<StringVector>(), "input files (ref and alt read counts)")
		("threshold,t", po::value<double>()->default_value(0.0), "threshold value for each strain")
		("lnorm,l", po::value<int>()->default_value(1), "the norm for the objective function")
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

    int nrThreads = vm["threads"].as<int>();
		int lnorm = vm["lnorm"].as<int>();
		double threshold = vm["threshold"].as<double>();
		std::cout << "threshold is " << threshold << std::endl;
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
		
    std::ofstream outFilteredInput(outputPrefix + "_filtered.tsv");
    outFilteredInput << filteredInput;
    outFilteredInput.close();
    
    std::ofstream outLog(outputPrefix + "_log.tsv");
		//outLog << "loc\tk\tLB\tUB" << std::endl;
		
		BoolMatrix Bmat = filteredInput.blowupBmat();
		
		std::ofstream outBmat(outputPrefix + "_blowupB.txt");
		for (int i = 0; i < Bmat.size(); ++i)
		{
			for (int j = 0; j < Bmat[i].size(); ++j)
			{
				outBmat << Bmat[i][j] << "\t";
			}
			outBmat << "\n";
		}
		
		int nrStrains = Bmat[0].size();
		
		Solver* psolve = nullptr;
		
		if (lnorm == 1)
		{
			psolve = new SolverRelax(filteredInput, Bmat, nrStrains, nrThreads, true, threshold);
		}
		else if (lnorm == 2)
		{
			psolve = new SolverRelax(filteredInput, Bmat, nrStrains, nrThreads, false, threshold);
		}
		
		if (psolve->solve())
		{
			Solution sol(psolve->getB(), psolve->getU(), filteredInput.getMutationDetails());
			
			std::ofstream outF(outputPrefix + "_F.txt");
			sol.writeSolF(filteredInput, outF);
			outF.close();
			
			std::ofstream outU(outputPrefix + "_U.txt");
			sol.writeSolU(filteredInput, outU);
			outU.close();
			
			std::ofstream outB(outputPrefix + "_B.txt");
			sol.writeSolB(filteredInput, outB);
			outB.close();
		}
		else
		{
			std::cout << "could not find solution!\n";
		}
		
  }
  catch (const std::runtime_error& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }

  return 0;
}
