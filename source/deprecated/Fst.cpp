//=======================================================================
// Author: Donovan Parks
//
// Copyright 2012 Donovan Parks
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

#include "Fst.hpp"

void Fst::Run(SplitSystem& splitSystem, Matrix& dissMatrix)
{
	dissMatrix.clear();

	// get outgroup sample
	uint outgroupId = INT_MAX;
	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(splitSystem.GetSampleName(i) == "Outgroup" || splitSystem.GetSampleName(i) == "outgroup")
			outgroupId = i;
	}

	// calculate within sample genetic distances
	std::vector<double> withinDiss;
	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(i == outgroupId)
			continue;

		std::vector<double> counts;
		double totalNumSeq;
		splitSystem.GetData(i, counts, totalNumSeq);

		double pairs = 0;
		double diss = 0;
		for(uint seqId1 = 0; seqId1 < counts.size(); ++seqId1)
		{
			if(counts.at(seqId1) == 0)
				continue;

			pairs += 0.5*counts.at(seqId1)*(counts.at(seqId1)-1);

			std::string seq1 = splitSystem.GetSequence(seqId1);

			for(uint seqId2 = 0; seqId2 < seqId1; ++seqId2)
			{
				if(counts.at(seqId2) == 0)
					continue;

				double weight = counts.at(seqId1)*counts.at(seqId2);
				pairs += weight;

				std::string seq2 = splitSystem.GetSequence(seqId2);

				diss += weight*PairwiseDifference(seq1, seq2);
			}
		}

		withinDiss.push_back(diss / pairs);
	}

	// calculate between sample genetic distances
	std::vector< std::vector<double> > betweenDiss;
	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(i == outgroupId)
			continue;

		std::vector<double> counts1;
		double totalNumSeq1;
		splitSystem.GetData(i, counts1, totalNumSeq1);

		std::vector<double> row;
		for(uint j = 0; j < i; ++j)
		{
			if(j == outgroupId)
				continue;

			std::vector<double> counts2;
			double totalNumSeq2;
			splitSystem.GetData(j, counts2, totalNumSeq2);

			double pairs = 0;
			double diss = 0;
			for(uint seqId1 = 0; seqId1 < counts1.size(); ++seqId1)
			{
				if(counts1.at(seqId1) == 0)
					continue;

				std::string seq1 = splitSystem.GetSequence(seqId1);

				for(uint seqId2 = 0; seqId2 < counts2.size(); ++seqId2)
				{
					if(counts2.at(seqId2) == 0)
						continue;

					double weight = counts1.at(seqId1)*counts2.at(seqId2);
					pairs += weight;

					std::string seq2 = splitSystem.GetSequence(seqId2);

					diss += weight*PairwiseDifference(seq1, seq2);
				}
			}

			row.push_back(diss / pairs);
		}

		betweenDiss.push_back(row);
	}

	// calculate Fst for each pair of samples
	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(i == outgroupId)
			continue;

		std::vector<double> row;
		for(uint j = 0; j < i; ++j)
		{
			if(j == outgroupId)
				continue;

			double within = 0.5 * (withinDiss.at(i) + withinDiss.at(j));
			double between = betweenDiss.at(i).at(j);

			row.push_back( (between - within) / between);
		}

		dissMatrix.push_back(row);
	}
}

uint Fst::PairwiseDifference(const std::string& seq1, const std::string& seq2)
{
	uint diff = 0;
	for(uint i = 0; i < seq1.size(); ++i)
	{
		if(seq1[i] != seq2[i])
			diff++;
	}

	return diff;
}