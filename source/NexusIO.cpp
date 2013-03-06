//=======================================================================
// Author: Donovan Parks
//
// Copyright 2009 Donovan Parks
//
// This file is part of Chameleon.
//
// Chameleon is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chameleon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Chameleon. If not, see <http://www.gnu.org/licenses/>.
//=======================================================================

#include "Precompiled.hpp"

#include "NexusIO.hpp"
#include "SplitSystem.hpp"
#include "Utils.hpp"


std::vector<std::string> &splitStr(const std::string &s, char delim, std::vector<std::string> &elems) 
{
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim))
		elems.push_back(item);

	return elems;
}


std::vector<std::string> splitStr(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	return splitStr(s, delim, elems);
}

bool NexusIO::Read(SplitSystem *const splitSystem, const std::string& filename)
{
	std::string nexusFile = filename;
	std::replace(nexusFile.begin(), nexusFile.end(), '\\', '/');

	std::ifstream textStream(nexusFile.c_str());
	if(!textStream.is_open())
		return false;

	// Read entire input stream
	// Note: Assumes TAXA block comes before the SPLITS block
	bool bReadTaxaBlock = false;
	std::string line;
	while(textStream.good())
  {		
		getline(textStream, line);
		if(line[0] != 'B')
			continue;

		if(line.find("BEGIN Taxa") != std::string::npos)
		{
			if(!ReadTaxaBlock(splitSystem, textStream))
				return false;

			bReadTaxaBlock = true;
		}
		else if(line.find("BEGIN Characters") != std::string::npos)
		{
			if(!ReadCharactersBlock(splitSystem, textStream))
				return false;
		}
		else if(line.find("BEGIN Trees") != std::string::npos)
		{
			if(!ReadTreesBlock(splitSystem, textStream))
				return false;
		}
		else if(line.find("BEGIN Splits")  != std::string::npos)
		{
			if(!bReadTaxaBlock)
			{
				std::cerr << "Unable to parse file. TAXA block is expected before SPLITS block.";
				return false;
			}

			if(!ReadSplitsBlock(splitSystem, textStream))
				return false;

			break;
		}
	}

	textStream.close();

	return true;
}

bool NexusIO::ReadTaxaBlock(SplitSystem *const splitSystem, std::ifstream& textStream)
{
	uint nexusId = 0;
	std::string line = "";
	std::set<std::string> seqNames;
	std::vector<std::string> missingSeqsInSampleFile;
	while(line.find("END; [Taxa]") == std::string::npos)
	{
		getline(textStream, line);

		if(line[0] == '[')
		{			
			int openingQuote = line.find('\'');
			int closingQuote = line.rfind('\'');
			std::string seqName = line.substr(openingQuote + 1, closingQuote - openingQuote - 1);

			// map sequence names to ids
			m_nexusIdToName[nexusId] = seqName;
			seqNames.insert(seqName);

			uint seqId;
			if(!splitSystem->GetSeqId(seqName, seqId))
				missingSeqsInSampleFile.push_back(seqName);

			nexusId++;
		}
	}

	// report missing sequences in sample file
	if(missingSeqsInSampleFile.size() > 0)
	{
		std::cout << "(Warning) The following taxa are in your split system, but not your sample file:" << std::endl;

		std::vector<std::string>::iterator iter;
		for(iter = missingSeqsInSampleFile.begin(); iter != missingSeqsInSampleFile.end(); ++iter)
			std::cout << *iter << std::endl;
	}

	// check for missing sequences in phylogeny
	splitSystem->CheckForMissingSeqs(seqNames);

	return true;
}

bool NexusIO::ReadCharactersBlock(SplitSystem *const splitSystem, std::ifstream& textStream)
{
	// move to sequence data
	std::string line = "";
	while(line.find("MATRIX") == std::string::npos)
		getline(textStream, line);

	// read sequence data
	while(line.find("END;") == std::string::npos)
		getline(textStream, line);

	return true;
}

bool NexusIO::ReadTreesBlock(SplitSystem *const splitSystem, std::ifstream& textStream)
{
	// skip tree block
	std::string line = "";
	while(line.find("END;") == std::string::npos)
		getline(textStream, line);

	return true;
}

bool NexusIO::ReadSplitsBlock(SplitSystem *const splitSystem, std::ifstream& textStream)
{	
	// determine number of splits
	std::string line = "";
	while(line.find("nsplits=") == std::string::npos)
		getline(textStream, line);

	uint indexOfSplitsCount = line.rfind('=');
	uint numSplits = atoi(line.substr(indexOfSplitsCount+1, line.rfind(';') - indexOfSplitsCount -1).c_str());

	// get to start of splits
	while(line.find("MATRIX") == std::string::npos)
		getline(textStream, line);

	// get list of outgroup sequences
	std::set<std::string> outgroupSeqs = splitSystem->GetOutgroupSeqs();

	// Sum of all sequences ids. Sequences ids go from 0 to one minus the number of sequences.
	uint numSeqs = splitSystem->GetNumIngroupSeqs();

	// determine all splits
	bool bError = false;
	for(int splitId = 0; splitId < (int)numSplits; ++splitId)
	{
		std::string splitLine;
		getline(textStream, splitLine);

		std::vector<std::string> rowData;
		splitStr(splitLine, '\t', rowData);

		double weight;
		std::string bipartitionStr;
		if(rowData.size() == 4)				// MedianNetwork
		{
			weight = fast_atof(rowData[2].c_str());
			bipartitionStr = rowData[3].substr(2, rowData[3].size()-3);
		}
		else if(rowData.size() == 3)	// NJ, NNet
		{
			weight = fast_atof(rowData[1].c_str());
			bipartitionStr = rowData[2].substr(2, rowData[2].size()-3);
		}
		else
		{
			std::cerr << "Unable to parse Nexus file. Invalid data in 'Splits' section";
			return false;
		}

		// bit array indicates the id of sequences on the left (1) and right (0) of the split
		std::vector<bool> split(numSeqs, false); 

		std::vector<std::string> ids;
		splitStr(bipartitionStr, ' ', ids);

		uint outgroupSeqOnLeft = 0;
		uint numSeqsOnLeft = 0;
		for(uint i = 0; i < ids.size(); ++i)
		{
			uint nexusId = atoi(ids.at(i).c_str())-1;
			std::string seqName = m_nexusIdToName[nexusId];

			if(outgroupSeqs.count(seqName) != 0)
				outgroupSeqOnLeft++;
			else
			{
				uint seqId;
				if(splitSystem->GetSeqId(seqName, seqId))
				{
					split[seqId] = true;
					numSeqsOnLeft++;		
				}
			}

		}

		uint outgroupSeqOnRight = outgroupSeqs.size() - outgroupSeqOnLeft;
		bool bOutgroupSeqOnLeft = (outgroupSeqOnLeft > 0);
		bool bOutgroupSeqOnRight = (outgroupSeqOnRight > 0);

		if(bOutgroupSeqOnLeft && bOutgroupSeqOnRight)
			continue;	// ignore splits with outgroup sequences in both subsets

		if(bOutgroupSeqOnLeft)
		{
			// swap split... always want outgroup seqs on the right
			for(uint i = 0; i < split.size(); ++i)
				split[i] = !split[i];

			numSeqsOnLeft = numSeqs - numSeqsOnLeft;
		}

		if(numSeqsOnLeft == numSeqs || numSeqsOnLeft == 0)
			continue;	// ignore 'outgroup' splits, or splits only containing taxa not found in the sample file
		
		splitSystem->AddSplit(Split(splitId, weight, split, false, true, numSeqsOnLeft, numSeqs-numSeqsOnLeft));
	}

	return !bError;
}
