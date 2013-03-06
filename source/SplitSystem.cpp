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

#include "SplitSystem.hpp"

#include "NexusIO.hpp"
#include "NewickIO.hpp"
#include "SampleIO.hpp"

bool SplitSystem::LoadData(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile, bool bVerbose)
{
	// read sample file
	if(!m_sampleIO.Read(sampleFile))
	{
		std::cerr << "Failed to read sample file: " << sampleFile;
		return false;
	}

	// read nexus file
	if(!nexusFile.empty())
	{
		NexusIO nexusIO;
		if(!nexusIO.Read(this, nexusFile))
		{
			std::cerr << "Failed to read Nexus file: " << nexusFile;
			return false;
		}
	}
	else if(!newickFile.empty())
	{
		NewickIO newickIO;
		if(!newickIO.Read(this, newickFile))
		{
			std::cerr << "Failed to read Newick file: " << nexusFile;
			return false;
		}
	}

	if(bVerbose)
	{
		std::cout << "  Number of samples: " << m_sampleIO.GetNumSamples() << std::endl;
		std::cout << "  Number of splits: " << GetNumSplits() << std::endl;
		std::cout << "  Total number of sequences: " << m_sampleIO.GetNumSeqs() << std::endl; 
		std::cout << "  Total number of ingroup sequences: " << m_sampleIO.GetNumIngroupSeqs() << std::endl; 
		std::cout << std::endl;
	}

	return true;
}

std::vector<double> SplitSystem::GetSampleData(uint sampleId, DATA_TYPE dataType)
{
	std::vector<double> data;
	data.resize(m_splits.size());

	std::vector<double> seqCount;
	double totalNumSeq;
	m_sampleIO.GetData(sampleId, seqCount, totalNumSeq);

	for(uint splitId = 0; splitId < m_splits.size(); ++splitId)
	{
		std::vector<uint> leftSeqIds = m_splits.at(splitId).GetLeftSequenceIds();
	
		for(uint i = 0; i < leftSeqIds.size(); ++i)
		{
			uint seqId = leftSeqIds.at(i);

			if(dataType == WEIGHTED_DATA)
				data.at(splitId) += seqCount.at(seqId) / totalNumSeq;
			else if(dataType == COUNT_DATA)
				data.at(splitId) += seqCount.at(seqId);
			else if(dataType == UNWEIGHTED_DATA)
			{
				if(seqCount.at(seqId) > 0)
					data.at(splitId) = 1;
			}
		}
	}

	return data;
}

void SplitSystem::CreateFromTree(Tree<Node>& tree)
{
	std::set<std::string> seqsToRemove = m_sampleIO.GetOutgroupSeqs();
	
	std::vector<std::string> seqsMissingInSampleFile;
	std::vector<Node*> leaves = tree.GetLeaves(tree.GetRootNode());
	for(uint i = 0; i < leaves.size(); ++i)
	{
		uint seqId;
		std::string name = leaves.at(i)->GetName();
		if(!m_sampleIO.GetSeqId(name, seqId))
		{
			seqsToRemove.insert(name);

			if(m_sampleIO.GetOutgroupSeqs().count(name) == 0)
				seqsMissingInSampleFile.push_back(name);
		}
	}

	// report missing sequences
	if(seqsMissingInSampleFile.size() > 0)
	{
		std::cout << "(Warning) The following sequences are in your tree file, but not your sample file:" << std::endl;

		std::vector<std::string>::iterator iter;
		for(iter = seqsMissingInSampleFile.begin(); iter != seqsMissingInSampleFile.end(); ++iter)
			std::cout << *iter << std::endl;
	}

	// project tree to ingroup sequences
	if(seqsToRemove.size() > 0)
		tree.Project(seqsToRemove);

	std::vector<Node*> nodes = tree.PostOrder(tree.GetRootNode());

	uint totalTaxa = m_sampleIO.GetNumIngroupSeqs();

	for(uint i = 0; i < nodes.size(); ++i)
	{
		Node* curNode = nodes.at(i);

		if(curNode->IsRoot())
			continue;

		// bit array indicates the id of sequences on the left (1) and right (0) of the split
		std::vector<bool> split(totalTaxa, false); 

		std::vector<Node*> leaves = curNode->GetLeaves();
		for(uint j = 0; j < leaves.size(); ++j)
		{
			std::string name = leaves.at(j)->GetName();
			uint id;
			if(m_sampleIO.GetSeqId(name, id))
				split[id] = true;
			else
			{
				assert(false);
				std::cerr << "(Bug) Unknown sequence name. Please report this bug." << std::endl;
			}
		}

		double weight = curNode->GetDistanceToParent();

		AddSplit(Split(i, weight, split, false, true, leaves.size(), totalTaxa-leaves.size()));
	}
}