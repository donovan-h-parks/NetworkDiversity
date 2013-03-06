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

#include "UnweightedUniFrac.hpp"

bool UnweightedUniFrac::Run(const SplitSystem& splitSystem, Matrix& dissMatrix)
{
	dissMatrix.clear();

	// get outgroup sample
	uint outgroupId = INT_MAX;
	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(splitSystem.GetSampleName(i) == "Outgroup" || splitSystem.GetSampleName(i) == "outgroup")
			outgroupId = i;
	}

	if(outgroupId == INT_MAX)
	{
		std::cout << "No outgroup defined." << std::endl;
		return false;
	}

	for(uint i = 0; i < splitSystem.GetNumSamples(); ++i)
	{
		if(i == outgroupId)
			continue;

		std::vector<double> row;
		for(uint j = 0; j < i; ++j)
		{
			if(j == outgroupId)
				continue;

			double shared = 0;
			double unique = 0;
			double root = 0;
			double external = 0;
			for(uint splitId = 0; splitId < splitSystem.GetNumSplits(); ++splitId)
			{
				double leftSeqCountI = splitSystem.GetCount(i, splitId);
				double rightSeqCountI = splitSystem.GetTotalNumSeq(i) - leftSeqCountI;

				double leftSeqCountJ = splitSystem.GetCount(j, splitId);				
				double rightSeqCountJ = splitSystem.GetTotalNumSeq(j) - leftSeqCountJ;

				double leftSeqCountOutgroup = splitSystem.GetCount(outgroupId, splitId);
				double rightSeqCountOutgroup = splitSystem.GetTotalNumSeq(outgroupId) - splitSystem.GetCount(outgroupId, splitId);
				
				if(leftSeqCountOutgroup > 0 && rightSeqCountOutgroup > 0)
				{
					std::cout << "problem" << std::endl;
					continue;
				}

				if(leftSeqCountOutgroup > 0)
				{
					double ingroupTaxaOnLeft = splitSystem.GetSplit(splitId).GetSizeLeftBipartition() - leftSeqCountOutgroup;

					if(rightSeqCountI > 0 && rightSeqCountJ == 0)
						unique += splitSystem.GetSplit(splitId).GetWeight();
					else if(rightSeqCountI == 0 && rightSeqCountJ > 0)
						unique += splitSystem.GetSplit(splitId).GetWeight();
					else if(rightSeqCountI > 0 && rightSeqCountJ > 0 && (leftSeqCountI+leftSeqCountJ) > 0)
						shared += splitSystem.GetSplit(splitId).GetWeight();
					else if(leftSeqCountI == 0 && leftSeqCountJ == 0 && ingroupTaxaOnLeft > 0)
						root += splitSystem.GetSplit(splitId).GetWeight();
					else if(rightSeqCountI == 0 && rightSeqCountJ == 0)
						external += splitSystem.GetSplit(splitId).GetWeight();
				}
				else if(rightSeqCountOutgroup > 0)
				{
					double ingroupTaxaOnRight = splitSystem.GetSplit(splitId).GetSizeRightBipartition() - rightSeqCountOutgroup;

					if(leftSeqCountI > 0 && leftSeqCountJ == 0)
						unique += splitSystem.GetSplit(splitId).GetWeight();
					else if(leftSeqCountI == 0 && leftSeqCountJ > 0)
						unique += splitSystem.GetSplit(splitId).GetWeight();
					else if(leftSeqCountI > 0 && leftSeqCountJ > 0 && (rightSeqCountI+rightSeqCountJ) > 0)
						shared += splitSystem.GetSplit(splitId).GetWeight();
					else if(rightSeqCountI == 0 && rightSeqCountJ == 0 && ingroupTaxaOnRight > 0)
						root += splitSystem.GetSplit(splitId).GetWeight();
					else if(leftSeqCountI == 0 && leftSeqCountJ == 0)
						external += splitSystem.GetSplit(splitId).GetWeight();
				}
			}

			row.push_back(unique / (unique + shared + root));
		}

		dissMatrix.push_back(row);
	}

	return true;
}