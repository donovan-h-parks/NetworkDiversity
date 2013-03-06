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

#include "Precompiled.hpp"

#include "DiversityCalculator.hpp"

bool DiversityCalculator::m_bWeighted;
std::vector< std::vector<double> > DiversityCalculator::m_dataVecRows;
std::vector< std::vector<double> > DiversityCalculator::m_dataVecCols;
std::vector<double> DiversityCalculator::m_minExtent;
std::vector<double> DiversityCalculator::m_maxExtent;
std::vector<double> DiversityCalculator::m_colSum;
std::vector<double> DiversityCalculator::m_rowLeafSum;
std::vector<double> DiversityCalculator::m_rowLeafSumSqrd;
std::vector<double> DiversityCalculator::m_weightedRowSum;
std::vector<double> DiversityCalculator::m_splitWeights;
double DiversityCalculator::m_totalSplitWeight;

DiversityCalculator::DiversityCalculator(SplitSystem& splitSystem, const std::string& calcStr, 
																								bool bWeighted, bool bCount, uint maxDataVecs, bool bVerbose)
	: m_maxDataVecs(maxDataVecs), m_bCount(bCount), m_bVerbose(bVerbose), m_bPhylogenetic(false), m_bGood(true), m_splitSystem(splitSystem)
{
	if(calcStr == "")
	{
		m_bGood = false;
		return;
	}

	std::clock_t divCalcStart = std::clock();

	m_bWeighted = bWeighted;

	if(!SetCalculator(calcStr))
	{
		m_bGood = false;
		return;
	}

	std::clock_t divCalcEnd = std::clock();

	if(m_bVerbose)
	{
		std::cout << "  Total time to initialize diversity calculator: " << ( divCalcEnd - divCalcStart ) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << std::endl;
	}
}

DiversityCalculator::~DiversityCalculator()
{

}

bool DiversityCalculator::SetCalculator(const std::string& calcStr)
{
	using namespace std::tr1::placeholders;

	bool bNeedColumnExtents = false;
	bool bNeedColumnSums = false;
	bool bNeedWeightedRowSums = false;
	bool bNeedTotalBranchLen = false;

	if(calcStr == "Bray-Curtis" || calcStr == "BC" || calcStr == "BrayCurtis")
		m_calculator = std::tr1::bind(&DiversityCalculator::BrayCurtis, _1, _2, _3, _4);
	else if(calcStr == "Canberra")
		m_calculator = std::tr1::bind(&DiversityCalculator::Canberra, _1, _2, _3, _4);
	else if(calcStr == "Coefficient of similarity" || calcStr == "CS" || calcStr == "CoefficientOfSimilarity")
		m_calculator = std::tr1::bind(&DiversityCalculator::CoefficientOfSimilarity, _1, _2, _3, _4);
	else if(calcStr == "Complete tree" || calcStr == "CT" || calcStr == "CompleteTree" || calcStr == "Complete Tree")
	{
		bNeedColumnExtents = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::CompleteTree, _1, _2, _3, _4);
	}
	else if(calcStr == "Euclidean")
		m_calculator = std::tr1::bind(&DiversityCalculator::Euclidean, _1, _2, _3, _4);
	else if(calcStr == "Gower")
	{
		bNeedColumnExtents = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::Gower, _1, _2, _3, _4);
	}
	else if(calcStr == "Kulczynski")
	{
		bNeedWeightedRowSums = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::Kulczynski, _1, _2, _3, _4);
	}
	else if(calcStr == "Lennon compositional difference" || calcStr == "Lennon" || calcStr == "LCD")
		m_calculator = std::tr1::bind(&DiversityCalculator::LennonCD, _1, _2, _3, _4);
	else if(calcStr == "Manhattan")
		m_calculator = std::tr1::bind(&DiversityCalculator::Manhattan, _1, _2, _3, _4);
	else if(calcStr == "Morisita-Horn"|| calcStr == "MH" || calcStr == "MorisitaHorn")
	{
		bNeedWeightedRowSums = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::MorisitaHorn, _1, _2, _3, _4);
	}
	else if(calcStr == "Soergel" || calcStr == "Ruzicka")
		m_calculator = std::tr1::bind(&DiversityCalculator::Soergel, _1, _2, _3, _4);
	else if(calcStr == "Tamas coefficient" || calcStr == "TC" || calcStr == "TamasCoefficient")
	{
		bNeedColumnExtents = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::TamasCoefficient, _1, _2, _3, _4);
	}
	else if(calcStr == "Weighted correlation" || calcStr == "WC" || calcStr == "WeightedCorrelation")
	{
		bNeedTotalBranchLen = true;
		bNeedWeightedRowSums = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::WeightedCorrelation, _1, _2, _3, _4);
	}
	else if(calcStr == "Yue-Clayton" || calcStr == "YC" || calcStr == "YueClayton")
		m_calculator = std::tr1::bind(&DiversityCalculator::YueClayton, _1, _2, _3, _4);
	else if(calcStr == "Sum")
		m_calculator = std::tr1::bind(&DiversityCalculator::Sum, _1, _2, _3, _4);
	else if(calcStr == "Extents")
	{
		bNeedColumnExtents = true;
		m_calculator = std::tr1::bind(&DiversityCalculator::Extents, _1, _2, _3, _4);
	}
	else
	{
		std::cerr << "Unknown calculator specified: " << calcStr << std::endl;
		return false;
	}

	// required to calculate intermediate terms
	GetSplitWeights();

	if(bNeedColumnExtents)
		CalculateColumnExtents();

	if(bNeedColumnSums)
		CalculateColumnSums();

	if(bNeedWeightedRowSums)
		CalculateWeightedRowSums();

	if(bNeedTotalBranchLen)
	{
		m_totalSplitWeight = 0;
		for(uint n = 0; n < m_splitWeights.size(); ++n)
			m_totalSplitWeight += m_splitWeights[n];
	}

	return true;
}

void DiversityCalculator::CalculateDataVectors(uint startIndex, uint numSamples, std::vector< std::vector<double> >& dataVec)
{
	std::clock_t startDataVecs = std::clock();
	
	// calculate data vector for each sample
	dataVec.clear();
	dataVec.reserve(numSamples);

	uint endIndex;
	if(m_splitSystem.IsOutgroup())
		endIndex = std::min<uint>(m_splitSystem.GetNumSamples()+1, startIndex+numSamples);
	else
		endIndex = std::min<uint>(m_splitSystem.GetNumSamples(), startIndex+numSamples);

	for(uint sampleId = startIndex; sampleId < endIndex; ++sampleId)
	{
		if(sampleId == m_splitSystem.GetOutgroupSampleId())
			continue;

		if(m_bWeighted && !m_bCount)
			dataVec.push_back(m_splitSystem.GetSampleData(sampleId, SplitSystem::WEIGHTED_DATA));
		else if(m_bWeighted && m_bCount)
			dataVec.push_back(m_splitSystem.GetSampleData(sampleId, SplitSystem::COUNT_DATA));
		else
			dataVec.push_back(m_splitSystem.GetSampleData(sampleId, SplitSystem::UNWEIGHTED_DATA));
	}

	std::clock_t endDataVecs = std::clock();
}

void DiversityCalculator::CalculateColumnExtents()
{
	std::clock_t extentsStart = std::clock();

	m_minExtent.clear();
	m_maxExtent.clear();
	m_minExtent.resize(m_splitSystem.GetNumSplits(), std::numeric_limits<double>::max());
	m_maxExtent.resize(m_splitSystem.GetNumSplits(), 0);

	for(uint i = 0; i < m_splitSystem.GetNumSamples(); i += m_maxDataVecs)
	{
		std::vector< std::vector<double> > data;
		CalculateDataVectors(i*m_maxDataVecs, m_maxDataVecs, data);
	
		for(uint j = 0; j < data.size(); ++j)
		{
			for(uint k = 0; k < data[j].size(); ++k)
			{
				if(data[j][k] < m_minExtent[k])
					m_minExtent[k] = data[j][k];

				if(data[j][k] > m_maxExtent[k])
					m_maxExtent[k] = data[j][k];
			}
		}
	}

	std::clock_t extentsEnd = std::clock();

	if(m_bVerbose)
	{
		std::cout << "  Time to calculate column extents: " << ( extentsEnd - extentsStart ) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << std::endl;
	}
}

void DiversityCalculator::CalculateColumnSums()
{
	std::clock_t colSumStart = std::clock();

	m_colSum.clear();
	m_colSum.resize(m_splitSystem.GetNumSplits(), 0);

	for(uint i = 0; i < m_splitSystem.GetNumSamples(); i += m_maxDataVecs)
	{
		std::vector< std::vector<double> > data;
		CalculateDataVectors(i*m_maxDataVecs, m_maxDataVecs, data);
	
		for(uint j = 0; j < data.size(); ++j)
		{
			for(uint k = 0; k < data[j].size(); ++k)
				m_colSum[k] = data[j][k];
		}
	}

	std::clock_t colSumEnd = std::clock();

	if(m_bVerbose)
	{
		std::cout << "  Time to calculate column sums: " << ( colSumEnd - colSumStart ) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << std::endl;
	}
}

void DiversityCalculator::CalculateWeightedRowSums()
{
	std::clock_t weightedRowSumStart = std::clock();

	m_weightedRowSum.clear();
	m_weightedRowSum.resize(m_splitSystem.GetNumSamples(), 0);

	for(uint i = 0; i < m_splitSystem.GetNumSamples(); i += m_maxDataVecs)
	{
		std::vector< std::vector<double> > data;
		CalculateDataVectors(i*m_maxDataVecs, m_maxDataVecs, data);
	
		for(uint j = 0; j < data.size(); ++j)
		{
			for(uint k = 0; k < data[j].size(); ++k)
				m_weightedRowSum[j] += m_splitWeights[k] * data[j][k];
		}
	}

	std::clock_t weightedRowSumEnd = std::clock();

	if(m_bVerbose)
	{
		std::cout << "  Time to calculate weighted row sums: " << ( weightedRowSumEnd - weightedRowSumStart ) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << std::endl;
	}
}

void DiversityCalculator::GetSplitWeights()
{
	m_splitWeights.clear();
	m_splitWeights.reserve(m_splitSystem.GetNumSplits());

	for(uint i = 0; i < m_splitSystem.GetNumSplits(); ++i)
		m_splitWeights.push_back(m_splitSystem.GetSplit(i).GetWeight());
}

bool DiversityCalculator::Dissimilarity(const std::string& dissFile)
{
	std::clock_t dissStart = std::clock();	

	// open dissimilarity file
	std::ofstream dissOut(dissFile.c_str());
	if(!dissOut.is_open())
	{
		std::cerr << "Unable to open dissimilarity matrix file: " << dissFile << std::endl;
		return false;
	}

	// get blocking information
	uint blockLen = m_maxDataVecs / 2;
	uint numBlocks = m_splitSystem.GetNumSamples() / blockLen;
	if(numBlocks*blockLen != m_splitSystem.GetNumSamples())
		++numBlocks;	// extra block if samples do not fit perfectly into blocks

	// calculate dissimilarity
	dissOut << m_splitSystem.GetNumSamples() << std::endl;

	double* partialDissMatrix = new double[blockLen*m_splitSystem.GetNumSamples()];

	double innerLoopTime = 0;
	for(uint row = 0; row < numBlocks; ++row)
	{
		CalculateDataVectors(row*blockLen, blockLen, m_dataVecRows);

		for(uint col = 0; col <= row; ++col)
		{
			CalculateDataVectors(col*blockLen, blockLen, m_dataVecCols);

			std::clock_t innerDissLoopStart = std::clock();	
			for(uint r = 0; r < m_dataVecRows.size(); ++r)
			{
				uint colStop = m_dataVecCols.size();
				if(row == 0)
					colStop = std::min<uint>(r, m_dataVecCols.size());

				for(uint c = 0; c < colStop; ++c)
				{
					double diss = m_calculator(m_dataVecRows[r], m_dataVecCols[c], r, c);
					partialDissMatrix[r*m_splitSystem.GetNumSamples() + col*blockLen + c] = diss;
				}
			}

			std::clock_t innerDissLoopEnd = std::clock();	
			innerLoopTime += (innerDissLoopEnd - innerDissLoopStart);
		}

		// write out partial dissimilarity matrix to file
		for(uint r = 0; r < m_dataVecRows.size(); ++r)
		{
			dissOut << m_splitSystem.GetSampleName(row*blockLen + r);

			for(uint c = 0; c < (row*blockLen + r); ++c)
				dissOut << '\t' << partialDissMatrix[r*m_splitSystem.GetNumSamples() + c];

			dissOut << std::endl;
		}
	}

	dissOut.close();

	delete[] partialDissMatrix;

	std::clock_t dissEnd = std::clock();

	if(m_bVerbose)
	{
		std::cout << std::endl;
		std::cout << "  Total time to calculate inner loop of dissimilarity matrix: " << innerLoopTime / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << "  Total time to calculate dissimilarity matrix: " << (dissEnd - dissStart) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << std::endl;
	}

	return true;
}

double DiversityCalculator::BrayCurtis(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double num = 0;
	double den = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		num += fabs(com1[n] - com2[n])*m_splitWeights[n];
		den += (com1[n] + com2[n])*m_splitWeights[n];
	}

	return num / den;
}

double DiversityCalculator::Canberra(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double diss = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		double den = com1[n] + com2[n];
		if(den != 0)
			diss += (fabs(com1[n] - com2[n]) / den)*m_splitWeights[n];
	}

	return diss;
}

double DiversityCalculator::CoefficientOfSimilarity(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double diss = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		double max = std::max<double>(com1[n],com2[n]);
		if(max > 0)
			diss += (fabs(com1[n]-com2[n])/max)*m_splitWeights[n];
	}

	return diss;
}

double DiversityCalculator::CompleteTree(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double num = 0;
	double den = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		num += fabs(com1[n] - com2[n])*m_splitWeights[n];
		den += (m_maxExtent[n] - m_minExtent[n])*m_splitWeights[n];
	}

	if(den == 0)
		return 1;

	return num / den;
}

double DiversityCalculator::Euclidean(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double diss = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		double d = com1[n] - com2[n];
		diss += m_splitWeights[n]*d*d;
	}

	return sqrt(diss);
}

double DiversityCalculator::Gower(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double diss = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		double d = m_maxExtent[n] - m_minExtent[n];
		if(d > 0)
			diss += (fabs(com1[n] - com2[n]) / (m_maxExtent[n] - m_minExtent[n]))*m_splitWeights[n];
	}

	return diss;
}

double DiversityCalculator::Kulczynski(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double sumMin = 0;
	for(uint n = 0; n < com1.size(); ++n)
		sumMin += std::min<double>(com1[n], com2[n])*m_splitWeights[n];

	return 1 - 0.5*(sumMin/m_weightedRowSum[i] + sumMin/m_weightedRowSum[j]);
}

double DiversityCalculator::LennonCD(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double A = 0;
	double B = 0;
	double C = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		A += std::min<double>(com1[n], com2[n])*m_splitWeights[n];
		B += (std::max<double>(com1[n], com2[n]) - com2[n])*m_splitWeights[n];
		C += (std::max<double>(com1[n], com2[n]) - com1[n])*m_splitWeights[n];
	}

	return std::min<double>(B, C) / (std::min<double>(B, C) + A);
}

double DiversityCalculator::Manhattan(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double diss = 0;
	for(uint n = 0; n < com1.size(); ++n)
		diss += fabs(com1[n] - com2[n])*m_splitWeights[n];

	return diss;
}

double DiversityCalculator::MorisitaHorn(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double prodSum = 0;
	double com1SumSqrd = 0;
	double com2SumSqrd = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		prodSum += com1[n]*com2[n]*m_splitWeights[n];

		com1SumSqrd += com1[n]*com1[n]*m_splitWeights[n];
		com2SumSqrd += com2[n]*com2[n]*m_splitWeights[n];
	}

	double num = 2*prodSum;
	double den = ((com1SumSqrd/(m_weightedRowSum[i]*m_weightedRowSum[i])) + (com2SumSqrd/(m_weightedRowSum[j]*m_weightedRowSum[j])))*m_weightedRowSum[i]*m_weightedRowSum[j];

	return 1.0 - num / den;
}

double DiversityCalculator::WeightedCorrelation(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double meanCol1 = m_weightedRowSum[i] / m_totalSplitWeight;
	double meanCol2 = m_weightedRowSum[j] / m_totalSplitWeight;

	double covXY = 0;
	double covX = 0;
	double covY = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		double diff1 = com1[n] - meanCol1;
		double diff2 = com2[n] - meanCol2;

		covXY += m_splitWeights[n]*diff1*diff2;
		covX += m_splitWeights[n]*diff1*diff1;
		covY += m_splitWeights[n]*diff2*diff2;
	}

	covXY /= m_totalSplitWeight;
	covX /= m_totalSplitWeight;
	covY /= m_totalSplitWeight;
	
	double denom = sqrt(covX*covY);
	if(denom == 0.0)
		return 0.0;

	double diss = 1.0 - covXY / denom;
	return diss;
}

double DiversityCalculator::Soergel(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double num = 0;
	double den = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		num += fabs(com1[n] - com2[n])*m_splitWeights[n];
		den += std::max<double>(com1[n], com2[n])*m_splitWeights[n];
	}

	return num / den;
}

double DiversityCalculator::TamasCoefficient(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double num = 0;
	double den = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		num += fabs(com1[n] - com2[n])*m_splitWeights[n];
		den += m_maxExtent[n]*m_splitWeights[n];
	}

	return num / den;
}

double DiversityCalculator::YueClayton(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double num = 0;
	double den = 0;
	for(uint n = 0; n < com1.size(); ++n)
	{
		num += com1[n]*com2[n]*m_splitWeights[n];

		double d = com1[n] - com2[n];
		den += (d*d + com1[n]*com2[n])*m_splitWeights[n];
	}

	return 1.0 - num / den;
}

double DiversityCalculator::Sum(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double sum = 0;
	for(uint n = 0; n < com1.size(); ++n)
		sum += (com1[n] + com2[n])*m_splitWeights[n];

	return sum;
}

double DiversityCalculator::Extents(const std::vector<double>& com1, const std::vector<double>& com2, uint i, uint j)
{
	double extents = 0;
	for(uint n = 0; n < com1.size(); ++n)
		extents += (m_maxExtent[n] - m_minExtent[n])*m_splitWeights[n];

	return extents;
}
