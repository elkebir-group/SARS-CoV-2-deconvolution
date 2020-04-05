/*
 * inputinstance.h
 *
 *  Created on: 3-apr-2020
 *      Author: M. El-Kebir
 */

#ifndef INPUTINSTANCE_H
#define INPUTINSTANCE_H

#include "utils.h"

/// Deconvolution input instance
class InputInstance
{
public:
  struct MutationDetails
  {
    /// Nucleotide position
    int _pos;
    /// Reference allele
    char _refAllele;
    /// altiant allele
    char _altAllele;
    /// Gene
    std::string _gene;
    /// Mutation type
    std::string _type;
    /// Amino acid substitution
    std::string _aminoAcidSub;
    /// Number of samples where mutation is subclonally present
    int _nrSubclonalSamples;
    /// Number of samples where mutation is clonally present
    int _nrClonalSamples;
    /// Number of samples where mutation is present
    int _nrSamples;
    /// Number of consensus sequence where mutation is present
    int _nrConsensusSamples;
  };
  
  enum MutationStatus
  {
    MutAbsent,
    MutClonal,
    MutSubclonal
  };
  
  /// Return mutation status
  /// @param i Mutation
  /// @param p Sample
  MutationStatus getMutationStatus(int i, int p) const
  {
    assert(0 <= i && i < _nrMutations);
    assert(0 <= p && p < _nrSamples);
    
    if (_vaf[i][p] >= 0.95 && _alt[i][p] >= 5)
    {
      return MutClonal;
    }
    else if (_vaf[i][p] >= 0.05 && _alt[i][p] >= 5)
    {
      return MutSubclonal;
    }
    else
    {
      return MutAbsent;
    }
  }
  
  /// Return number of samples
  int getNrSamples() const
  {
    return _nrSamples;
  }
  
  /// Return number of mutations
  int getNrMutations() const
  {
    return _nrMutations;
  }
  
  /// Return sample labels
  const StringVector& getSamples() const
  {
    return _samples;
  }
  
  /// Return altiant allele frequency
  /// @param i Mutation
  /// @param p Sample
  double getVaf(int i, int p) const
  {
    assert(0 <= i && i < _nrMutations);
    assert(0 <= p && p < _nrSamples);
    
    return _vaf[i][p];
  }
  
  /// Return alt read count
  /// @param i Mutation
  /// @param p Sample
  double getAltCount(int i, int p) const
  {
    assert(0 <= i && i < _nrMutations);
    assert(0 <= p && p < _nrSamples);
    
    return _alt[i][p];
  }
  
  /// Return reference read count
  /// @param i Mutation
  /// @param p Sample
  double getRefCount(int i, int p) const
  {
    assert(0 <= i && i < _nrMutations);
    assert(0 <= p && p < _nrSamples);
    
    return _ref[i][p];
  }
  
  /// Return mutation details
  /// @param i Mutation
  const MutationDetails& getMutationDetails(int i) const
  {
    assert(0 <= i && i < _nrMutations);
    return _mutDetails[i];
  }
  
  /// Read input files
  /// @param inRef Reference reads
  /// @param inAlt Alternate reads
  void read(std::istream& inRef, std::istream& inAlt);
  
  /// Filter out local mutations
  /// @param newToOld Map from new indices to old indices
  /// @param oldToNew Map from old indices to new indices (-1 indicates unmapped)
  InputInstance filterMutations(IntVector& newToOld,
                                     IntVector& oldToNew) const;
  
  /// Filter out samples with no mutations
  /// @param newToOld Map from new indices to old indices
  /// @param oldToNew Map from old indices to new indices (-1 indicates unmapped)
  InputInstance filterSamples(IntVector& newToOld, IntVector& oldToNew) const;
  
private:
  typedef std::vector<MutationDetails> MutationDetailsVector;
  
  /// Number of bulk samples
  int _nrSamples;
  /// Number of mutations
  int _nrMutations;
  /// Mutation details
  MutationDetailsVector _mutDetails;
  /// Sample names
  StringVector _samples;
  /// altiant allele frequencies
  DoubleMatrix _vaf;
  /// Reference read count
  IntMatrix _ref;
  /// altiant read count
  IntMatrix _alt;
  
  friend std::ostream& operator<<(std::ostream& out, const InputInstance& input);
  friend std::istream& operator>>(std::istream& in, InputInstance& input);
};

#endif // INPUTINSTANCE_H
