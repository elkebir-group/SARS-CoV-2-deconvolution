/*
 * utils.h
 *
 *  Created on: 19-oct-2017
 *      Author: M. El-Kebir
 */

#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <iostream>
#include <set>
#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

typedef std::map<std::string, int> StringToIntMap;
typedef std::vector<std::string> StringVector;
typedef std::list<std::string> StringList;
typedef std::set<std::string> StringSet;
typedef std::set<StringSet> StringSetSet;
typedef std::set<int> IntSet;
typedef std::list<int> IntList;
typedef std::vector<int> IntVector;
typedef std::vector<IntVector> IntMatrix;
typedef std::vector<IntMatrix> IntTensor;
typedef std::vector<bool> BoolVector;
typedef std::vector<BoolVector> BoolMatrix;
typedef std::vector<BoolMatrix> BoolTensor;
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleMatrix;
typedef std::vector<DoubleMatrix> DoubleTensor;
typedef std::pair<std::string, std::string> StringPair;
typedef std::list<StringPair> StringPairList;
typedef std::pair<int, double> IntDoublePair;
typedef std::vector<IntDoublePair> IntDoublePairVector;
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
typedef std::vector<IntSet> IntSetVector;
typedef std::pair<int, IntPair> IntTriple;

/// Get line from stream in a platform independent manner
std::istream& getline(std::istream& is, std::string& t);

/// Random number generator
extern std::mt19937 g_rng;

/// Current line number
extern int g_lineNumber;

/// Get current line number as a string
std::string getLineNumber();

#endif // UTILS_H
