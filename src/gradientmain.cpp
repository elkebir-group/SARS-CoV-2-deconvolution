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
  SolverGradient::Param param;
  
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("strains,k", po::value<int>(), "number of strains")
    ("lambda,l", po::value<double>()->default_value(param._lambda), "rate of increase of lambda")
	("lambdaInit,linit", po::value<double>()->default_value(param._lambdaInit), "initial value of lambda")
    ("initB,B", po::value<std::string>(), "genotype matrix for initialization")
    ("initU,U", po::value<std::string>(), "mixture matrix for initialization")
    ("maxIter,m", po::value<int>()->default_value(param._maxIter), "maximum number of iterations")
    ("eps,e", po::value<double>()->default_value(param._epsilon), "termination condition")
		("filter,f", po::bool_switch()->default_value(false), "filtering flag")
    ("seed,s", po::value<int>()->default_value(0), "random number generator seed")
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

    param._nrThreads = vm["threads"].as<int>();
    param._nrStrains = vm["strains"].as<int>();
    param._epsilon = vm["eps"].as<double>();
    param._maxIter = vm["maxIter"].as<int>();
    param._lambda = vm["lambda"].as<double>();
		param._lambdaInit = vm["lambdaInit"].as<double>();
    const bool filter = vm["filter"].as<bool>();

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
    
    int seed = vm["seed"].as<int>();
    g_rng = std::mt19937(seed);

		InputInstance filteredInput;
		if (!filter)
		{
			filteredInput = input;
		}
		else
		{
			std::cout << "performing filtering\n";
			filteredInput = input.filter();
		}

    std::ofstream outFilteredInput(outputPrefix + "_filtered.tsv");
    outFilteredInput << filteredInput;
    outFilteredInput.close();
    
    SolverGradient solver(filteredInput, param, outputPrefix);
    
    if (vm.count("initB"))
    {
      std::ifstream in(vm["initB"].as<std::string>().c_str());
      Solution sol;
      sol.readSolB(in);
      solver.initB(sol.getDoubleB());
    }
    
    if (vm.count("initU"))
    {
      std::ifstream in(vm["initU"].as<std::string>().c_str());
      Solution sol;
      sol.readSolU(in);
      solver.initU(sol.getU());
    }
    
    if (solver.solve())
		{
      Solution sol(solver.getDoubleB(), solver.getU(), filteredInput.getMutationDetails());
			
			std::ofstream outF(outputPrefix + "_F.txt");
			sol.writeSolF(filteredInput, outF);
			outF.close();
			
			std::ofstream outU(outputPrefix + "_U.txt");
			sol.writeSolU(filteredInput, outU);
			outU.close();
			
			std::ofstream outB(outputPrefix + "_B.txt");
			sol.writeSolB(filteredInput, outB);
			outB.close();
      
      std::ofstream outDoubleB(outputPrefix + "_doubleB.txt");
      sol.writeSolDoubleB(filteredInput, outDoubleB);
      outDoubleB.close();
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
  catch (const po::error& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }

  return 0;
}

