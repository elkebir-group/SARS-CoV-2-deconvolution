/*
 * solution.h
 *
 *  Created on: 8-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef SOLUTION_H
#define SOLUTION_H

#include "inputinstance.h"

class Solution
{
public:
  /// Default constructor
  Solution();
  
  /// Constructor
  /// @param doubleB Genotype matrix
  /// @param U Mixture matrix
  /// @param mutDetails Mutation details
  Solution(const DoubleMatrix& doubleB, const DoubleMatrix& U,
           const InputInstance::MutationDetailsVector& mutDetails);
  
  /// Constructor
  /// @param B Genotype matrix
  /// @param U Mixture matrix
  /// @param mutDetails Mutation details
  Solution(const BoolMatrix& B, const DoubleMatrix& U,
           const InputInstance::MutationDetailsVector& mutDetails);
  
  /// Compute L2 distance
  /// @param input Input instance
  double computeDistanceL2(const InputInstance& input) const;
  
  /// Compute Bayesian Information Criterion
  /// @param input Input instance
  double getBIC(const InputInstance& input) const;
  
  int getNrPresentMutations() const;
  
  /// Return genotype of specified strain
  /// @param j Strain
  BoolVector getGenotype(int j) const;
  
  /// Return extended genotype of specified strain
  /// @param input Input instance
  /// @param j Strain
  BoolVector getExtGenotype(const InputInstance& input,
                            int j) const;
  
  /// Return number of strains
  int getNrStrains() const
  {
    //assert(_U.size() > 0);
    assert(_B.size() > 0);
    //assert(_B[0].size() == _U.size());
    
    return _B[0].size();
  }
  
  /// Return number of samples
  int getNrSamples() const
  {
    assert(_U.size() > 0);
    return _U[0].size();
  }
  
  /// Return number of mutations
  int getNrMutations() const
  {
    return _B.size();
  }
  
  /// Return genotype matrix B
  const DoubleMatrix& getDoubleB() const
  {
    return _doubleB;
  }
  
  /// Return genotype matrix B
  const BoolMatrix& getB() const
  {
    return _B;
  }
  
  /// Return mixture matrix U
  const DoubleMatrix& getU() const
  {
    return _U;
  }
  
  /// Return frequency matrix F
  const DoubleMatrix& getF() const
  {
    return _F;
  }
  
  /// Write frequency matrix
  /// @param input Input instance
  /// @param out Output stream
  void writeSolF(const InputInstance& input,
                 std::ostream& out) const;
  
  /// Write mixture matrix
  /// @param input Input instance
  /// @param out Output stream
  void writeSolU(const InputInstance& input,
                 std::ostream& out) const;
  
  /// Write genotype matrix
  /// @param input Input instance
  /// @param out Output stream
  void writeSolB(const InputInstance& input,
                 std::ostream& out) const;
  
  /// Write genotype matrix (double)
  /// @param input Input instance
  /// @param out Output stream
  void writeSolDoubleB(const InputInstance& input,
                       std::ostream& out) const;
  
  /// Read solution
  /// @param inB Input stream for B
  /// @param inU Input stream for U
  void readSol(std::istream& inB,
               std::istream& inU);
  
  /// Read mixture matrix
  /// @param in Input stream
  void readSolU(std::istream& in);
  
  /// Read genotype matrix
  /// @param in Input stream
  void readSolB(std::istream& in);
  
private:
  /// Genotype matrix (mutations by strains)
  BoolMatrix _B;
  /// Genotype matrix as double (mutations by strains)
  DoubleMatrix _doubleB;
  /// Mixture matrix (strains by samples)
  DoubleMatrix _U;
  /// Inferred frequences (mutations by samples)
  DoubleMatrix _F;
  /// Mutation details
  InputInstance::MutationDetailsVector _mutDetails;
};

#endif // SOLUTION_H
