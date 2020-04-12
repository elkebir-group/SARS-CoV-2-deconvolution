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
  /// Default constructor
  InputInstance();
  
  /// Constructor
  /// @param F Frequency matrix (mutation by sample)
  /// @param depth Number of reads
  InputInstance(DoubleMatrix F,
                int depth);
  
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
  
  typedef std::vector<MutationDetails> MutationDetailsVector;
  
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
  
  /// Return sample location
  const StringVector& getSampleLocations() const
  {
    return _sampleLocations;
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
  
  /// Return mutation details
  const MutationDetailsVector getMutationDetails() const
  {
    return _mutDetails;
  }
  
  /// Read input files
  /// @param inRef Reference reads
  /// @param inAlt Alternate reads
  void read(std::istream& inRef, std::istream& inAlt);
  
  /// Write input files
  /// @param outRef Reference reads
  /// @param outAlt Alternate reads
  void write(std::ostream& outRef, std::ostream& outAlt) const;
  
  /// Filter samples and mutations
  InputInstance filter() const;
  
  /// Filter out local mutations
  /// @param newToOld Map from new indices to old indices
  /// @param oldToNew Map from old indices to new indices (-1 indicates unmapped)
  InputInstance filterMutations(IntVector& newToOld,
                                IntVector& oldToNew) const;
  
  /// Filter out samples with no mutations
  /// @param newToOld Map from new indices to old indices
  /// @param oldToNew Map from old indices to new indices (-1 indicates unmapped)
  InputInstance filterSamples(IntVector& newToOld, IntVector& oldToNew) const;
  
  /// Filter out samples by location
  /// @param location Location
  /// @param newToOld Map from new indices to old indices
  /// @param oldToNew Map from old indices to new indices (-1 indicates unmapped)
  InputInstance filterSamplesByLocation(const std::string& location,
                                        IntVector& newToOld,
                                        IntVector& oldToNew) const;
  
  InputInstance filterSamples(const IntSet& sampleIndices) const;
  
  typedef std::map<std::string, InputInstance> InputInstanceMap;
  
  InputInstanceMap splitSamplesByLocation() const;

	BoolMatrix blowupBmat();
	
protected:
  template <typename T>
  void writeTemp(const T& matrix, std::ostream& out) const;
  
  void initSampleLocations();
  
private:
  /// Number of bulk samples
  int _nrSamples;
  /// Number of mutations
  int _nrMutations;
  /// Mutation details
  MutationDetailsVector _mutDetails;
  /// Sample names
  StringVector _samples;
  /// Sample locations
  StringVector _sampleLocations;
  /// altiant allele frequencies
  DoubleMatrix _vaf;
  /// Reference read count
  IntMatrix _ref;
  /// altiant read count
  IntMatrix _alt;
	
  friend std::ostream& operator<<(std::ostream& out, const InputInstance& input);
  friend std::istream& operator>>(std::istream& in, InputInstance& input);
};

template <typename T>
void InputInstance::writeTemp(const T& matrix, std::ostream& out) const
{
  out << "pos" << "\t"
      << "ref" << "\t"
      << "alt" << "\t"
      << "gene" << "\t"
      << "N/S" << "\t"
      << "AA" << "\t"
      << "nSRAsubclonal" << "\t"
      << "nSRAclonal" << "\t"
      << "nSRA" << "\t"
      << "nConsensus";
  
  for (const std::string& sample : _samples)
  {
    out << "\t" << sample;
  }
  
  out << std::endl;
  
  const int n = _nrMutations;
  const int m = _nrSamples;
  for (int i = 0; i < n; ++i)
  {
    const InputInstance::MutationDetails& mutDetails_i = _mutDetails[i];
    out << mutDetails_i._pos << "\t"
        << mutDetails_i._refAllele << "\t"
        << mutDetails_i._altAllele << "\t"
        << mutDetails_i._gene << "\t"
        << mutDetails_i._type << "\t"
        << mutDetails_i._aminoAcidSub << "\t"
        << mutDetails_i._nrSubclonalSamples << "\t"
        << mutDetails_i._nrClonalSamples << "\t"
        << mutDetails_i._nrSamples << "\t"
        << mutDetails_i._nrConsensusSamples;
    
    for (int p = 0; p < m; ++p)
    {
      out << "\t";

      const double x_ip = matrix[i][p];
      if (std::isnan(x_ip))
      {
        out << "NULL";
      }
      else
      {
        out << x_ip;
      }
    }
    
    out << std::endl;
  }
}

#endif // INPUTINSTANCE_H
