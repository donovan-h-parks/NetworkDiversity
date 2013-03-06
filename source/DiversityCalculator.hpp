//=======================================================================
// Author: Donovan Parks
//
// Copyright 2011 Donovan Parks
//
// This file is part of ExpressBetaDiversity.
//
// ExpressBetaDiversity is free software: you can redistribute it 
// and/or modify it under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 3 of 
// the License, or (at your option) any later version.
//
// ExpressBetaDiversity is distributed in the hope that it will be 
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ExpressBetaDiversity. If not, see 
// <http://www.gnu.org/licenses/>.
//=======================================================================

#ifndef _DIVERSITY_CALCULATOR_
#define _DIVERSITY_CALCULATOR_

#include "Precompiled.hpp"

#include "SplitSystem.hpp"
/**
 * @brief Measure beta-diversity with a variety of calculators.
 */
class DiversityCalculator
{
public:		
	/** Constructor. */
	DiversityCalculator(SplitSystem& splitSystem, const std::string& calcStr, 
													bool bWeighted, bool bCount = false, uint maxDataVecs = 1000, bool bVerbose = false);

	/** Destructor. */
	~DiversityCalculator();

	/** Check good flag. */
	bool IsGood() const { return m_bGood; }

	/** Calculate dissimilarity between all pairs of samples. */
	bool Dissimilarity(const std::string& dissFile);

private:
	/** Set desired calculator. */
	bool SetCalculator(const std::string& calcStr);

	/** Calculate data vectors . */
	void CalculateDataVectors(uint startIndex, uint numSamples, std::vector< std::vector<double> >& dataVec);

	/** Calculate minimum and maximum value of each column in the sample data matrix. */
	void CalculateColumnExtents();

	/** Calculate sum of each column in the sample data matrix. */
	void CalculateColumnSums();

	/** Calculate branch length weighted sum for each row in the data matrix. */
	void CalculateWeightedRowSums();

	/** Get weight of each split. */
	void GetSplitWeights();

	static double BrayCurtis(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Canberra(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double CoefficientOfSimilarity(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double CompleteTree(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Euclidean(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Gower(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Kulczynski(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double LennonCD(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Manhattan(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double MorisitaHorn(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Soergel(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double TamasCoefficient(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double WeightedCorrelation(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double YueClayton(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);

	static double Sum(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	static double Extents(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j);
	
private:
	typedef std::tr1::function<double (const std::vector<double>&, const std::vector<double>&, uint, uint)> CalculatorFunc;

	/** Split system to calculate beta diversity over. */
	SplitSystem& m_splitSystem;

	/** Function object indicating calculator to use. */
	CalculatorFunc m_calculator;

	/** Flag indicating if all is good in the world. */
	bool m_bGood;

	/** Maximum number of data vectors to have in memory at once. */
	uint m_maxDataVecs;

	/** Flag indicating if weighted vectors are to be generated. */
	static bool m_bWeighted;

	/** Flag indicating if count data should be used or if it should be normalized to relative proportions. */
	bool m_bCount;

	/** Flag indicating if phylogenetic vectors are to be generated. */
	bool m_bPhylogenetic;

	/** Flag indicating if program execution information should be printed. */
	bool m_bVerbose;

	/** Weight associated with each split/column. */
	static std::vector<double> m_splitWeights;

	/** Total split weight in split system. */
	static double m_totalSplitWeight;

	/** Data vectors for current rows in dissimilarity matrix being processed. */ 
	static std::vector< std::vector<double> > m_dataVecRows;

	/** Data vectors for current columns in dissimilarity matrix being processed. */ 
	static std::vector< std::vector<double> > m_dataVecCols;

	/** Minimum value in each column of data matrix. */
	static std::vector<double> m_minExtent;

	/** Maximum value in each column of data matrix. */
	static std::vector<double> m_maxExtent;

	/** Sum of each column in the data matrix. */
	static std::vector<double> m_colSum;

	/** Sum of leaf node proportions along each row of the data matrix. */
	static std::vector<double> m_rowLeafSum;

	/** Sum of squarded leaf node proportions along each row of the data matrix. */
	static std::vector<double> m_rowLeafSumSqrd;

	/** Sum of each row weighted by branch length in the data matrix. */
	static std::vector<double> m_weightedRowSum;
};

#endif

