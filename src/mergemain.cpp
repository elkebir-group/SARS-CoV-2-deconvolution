/*
 * mergeemain.cpp
 *
 *  Created on: 8-apr-2020
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "inputinstance.h"
#include "solution.h"
#include "solutionset.h"
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

StringSet getSampleLocationSet(const InputInstance& input)
{
  StringSet res;
  
  for (const std::string& sample : input.getSampleLocations())
  {
    res.insert(sample);
  }
  
  return res;
}

int main(int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("strains,k", po::value<int>(), "maximum number of strains")
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
    
    SolutionSet solSet(input);
    solSet.populate(outputPrefix, vm["strains"].as<int>());
    
    Solution mergedSol = solSet.merge();
    
    std::ofstream outB(outputPrefix + "_merged_B.tsv");
    mergedSol.writeSolB(solSet.getInput(), outB);
    outB.close();
    
    std::ofstream outU(outputPrefix + "_merged_U.tsv");
    mergedSol.writeSolU(solSet.getInput(), outU);
    outU.close();
    
    std::ofstream outF(outputPrefix + "_merged_F.tsv");
    mergedSol.writeSolF(solSet.getInput(), outF);
    outF.close();
  }
  catch (const std::runtime_error& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  catch (const boost::program_options::error& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  
  return 0;
}
