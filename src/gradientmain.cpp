/*
 * gradientmain.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "inputinstance.h"
#include "solution.h"
#include <fstream>
#include <boost/program_options.hpp>
#include "solvergradient.h"

namespace po = boost::program_options;

int main(int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("strains,k", po::value<int>(), "number of strains")
    ("restarts,N", po::value<int>()->default_value(50), "number of restarts")
    ("eps,e", po::value<double>()->default_value(0.01), "termination condition")
    ("threads,T", po::value<int>()->default_value(1), "number of threads")
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

    const int nrThreads = vm["threads"].as<int>();
    const int nrStrains = vm["strains"].as<int>();
    const int nrRestarts = vm["restarts"].as<int>();
    const double epsilon = vm["eps"].as<double>();
    
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

    InputInstance filteredInput = input;//.filter();

    std::ofstream outFilteredInput(outputPrefix + "_filtered.tsv");
    outFilteredInput << filteredInput;
    outFilteredInput.close();
    
    SolverGradient solver(filteredInput, nrStrains, nrRestarts, nrThreads, epsilon);
    if (solver.solve())
		{
      Solution sol(solver.getB(), solver.getU(), filteredInput.getMutationDetails());
			
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

