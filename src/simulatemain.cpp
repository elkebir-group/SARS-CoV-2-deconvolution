/*
 * simulatemain.cpp
 *
 *  Created on: 11-apr-2020
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "inputinstance.h"
#include "solution.h"
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace po = boost::program_options;
namespace ub = boost::numeric::ublas;

int main(int argc, char** argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("strains,k", po::value<int>()->default_value(10), "number of strains")
    ("missing", po::value<double>()->default_value(0), "missing mutation rate")
    ("samples,m", po::value<int>()->default_value(50), "number of samples")
    ("expstrains,K", po::value<int>()->default_value(3), "number of expected strains per sample")
    ("expmutations,N", po::value<int>()->default_value(5), "number of expected mutations per strain")
    ("mutations,n", po::value<int>()->default_value(100), "number of mutations")
    ("depth,d", po::value<int>()->default_value(1000), "number of reads")
    ("seed,s", po::value<int>()->default_value(0), "random number generator seed")
    ("output,o", po::value<std::string>()->default_value("out"), "output prefix")
  ;
  
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
    
    const int nrStrains = vm["strains"].as<int>();
    const int nrExpStrains = vm["expstrains"].as<int>();
    const int nrSamples = vm["samples"].as<int>();
    const int nrMutations = vm["mutations"].as<int>();
    const int nrExpMutations = vm["expmutations"].as<int>();
    const int depth = vm["depth"].as<int>();
    const std::string outputPrefix = vm["output"].as<std::string>();
    const double missing = vm["missing"].as<double>();
    
    int seed = vm["seed"].as<int>();
    g_rng = std::mt19937(seed);
    
    // generate B
    std::uniform_real_distribution<> unif(0, 1);
    
    ub::matrix<double> B(nrMutations, nrStrains);
    for (int j = 0; j < nrStrains; ++j)
    {
      double prob = (double)nrExpMutations / nrMutations;
      for (int i = 0; i < nrMutations; ++i)
      {
        B(i, j) = unif(g_rng) < prob ? 1. : 0.;
      }
    }
    
    // generate U
    std::poisson_distribution<> pois(nrExpStrains - 1);
    
    IntVector strains;
    for (int j = 0; j < nrStrains; ++j)
    {
      strains.push_back(j);
    }
    
    // https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution
    // Symmetric Dirichlet with concentration parameter alpha = 1
    std::gamma_distribution<> gamma(1, 1);
    
    ub::matrix<double> U(nrStrains, nrSamples, 0);
    
    for (int p = 0; p < nrSamples; ++p)
    {
      // 1. decide on how many strains for this sample using a Poisson
      const int nrStrains_p = std::min(1 + pois(g_rng), nrStrains);
      
      // 2. pick strains
      std::random_shuffle(strains.begin(), strains.end());
      
      // 3. draw from Dirichlet
      double sum = 0;
      for (int jj = 0; jj < nrStrains_p; ++jj)
      {
        int j = strains[jj];
        U(j, p) = gamma(g_rng);
        sum += U(j, p);
      }
      for (int jj = 0; jj < nrStrains_p; ++jj)
      {
        int j = strains[jj];
        U(j, p) /= sum;
      }
    }
    
    ub::matrix<double> M(nrMutations, nrSamples, 0);
    for (int i = 0; i < nrMutations; ++i)
    {
      for (int p = 0; p < nrSamples; ++p)
      {
        if (unif(g_rng) < missing)
        {
          M(i, p) = -1;
        }
      }
    }
    
    // generate F
    ub::matrix<double> F = ub::prod(B, U);
    
    DoubleMatrix stdF(nrMutations, DoubleVector(nrSamples, 0));
    for (int i = 0; i < nrMutations; ++i)
    {
      for (int p = 0; p < nrSamples; ++p)
      {
        if (M(i,p) ==  -1)
        {
          stdF[i][p] = -1;
        }
        else
        {
          stdF[i][p] = F(i, p);
        }
      }
    }
    
    BoolMatrix stdB(nrMutations, BoolVector(nrStrains, false));
    for (int i = 0; i < nrMutations; ++i)
    {
      for (int j = 0; j < nrStrains; ++j)
      {
        stdB[i][j] = B(i, j) > 0.5;
      }
    }
    
    DoubleMatrix stdU(nrStrains, DoubleVector(nrSamples, 0.));
    for (int j = 0; j < nrStrains; ++j)
    {
      for (int p = 0; p < nrSamples; ++p)
      {
        stdU[j][p] = U(j, p);
      }
    }
    
    InputInstance input(stdF, depth);
    
    std::ofstream outAlt(outputPrefix + "_alt.tsv");
    std::ofstream outRef(outputPrefix + "_ref.tsv");
    input.write(outRef, outAlt);
    outRef.close();
    outAlt.close();
    
    std::ofstream outVaf(outputPrefix + "_vaf.tsv");
    outVaf << input;
    outVaf.close();
    
    Solution sol(stdB, stdU, input.getMutationDetails());
    
    std::ofstream outB(outputPrefix + "_B.tsv");
    sol.writeSolB(input, outB);
    outB.close();
    
    std::ofstream outU(outputPrefix + "_U.tsv");
    sol.writeSolU(input, outU);
    outU.close();
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
