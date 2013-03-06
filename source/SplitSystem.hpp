//=======================================================================
// Author: Donovan Parks
//
// Copyright 2009 Donovan Parks
//
// This file is part of Chameleon.

// Chameleon is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Chameleon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Chameleon. If not, see <http://www.gnu.org/licenses/>.
//=======================================================================

#ifndef _SPLIT_SYSTEM_
#define _SPLIT_SYSTEM_

#include "Precompiled.hpp"

#include "Split.hpp"
#include "SampleIO.hpp"

#include "Tree.hpp"
#include "Node.hpp"

/**
 * @class SplitSystem
 * @brief Holds data and functions for a collection of splits.
 */
class SplitSystem
{
public:
	enum DATA_TYPE { WEIGHTED_DATA, COUNT_DATA, UNWEIGHTED_DATA };

public:
	/** Constructor. */
	SplitSystem() {}

	/** Destructor. */
	~SplitSystem() {}

	/** Load input data. */
	bool LoadData(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile, bool bVerbose = false);

	/** Create split system from a tree. */
	void CreateFromTree(Tree<Node>& tree);

	/** Check for sequences in sample file not in the phylogeny */
	void CheckForMissingSeqs(const std::set<std::string>& seqsInPhylogeny) { m_sampleIO.CheckForMissingSeqs(seqsInPhylogeny); }

	/** Get number of splits. */
	uint GetNumSplits() const { return m_splits.size(); }

	/** 
	 * @brief Get a desired split.
	 * @param index Index of desired split.
	 * @return Split at given index.
	 */
	const Split& GetSplit(uint index) const { return m_splits.at(index); }

	/** Add split. */
	void AddSplit(const Split& split) { m_splits.push_back(split); }

	/** Get number of sequences (including any outgroup or missing sequences). */
	uint GetNumSeqs() const { return m_sampleIO.GetNumSeqs(); }

	/** Get number of ingroup sequences (excludes any outgroup or missing sequences). */
	uint GetNumIngroupSeqs() const { return m_sampleIO.GetNumIngroupSeqs(); }

	/** Get number of samples. */
	uint GetNumSamples() const { return m_sampleIO.GetNumSamples(); } 

	/** Get name of sample. */
	std::string GetSampleName(uint sampleId) const { return m_sampleIO.GetSampleName(sampleId); }

	/** Get sequence id.*/
	bool GetSeqId(const std::string& name, uint& seqId) { return m_sampleIO.GetSeqId(name, seqId); }

	/** Get data from specified sample. */
	std::vector<double> GetSampleData(uint sampleId, DATA_TYPE dataType);

	/** Check if there is an outgroup. */
	bool IsOutgroup() const { return m_sampleIO.IsOutgroup(); }

	/** Get outgroup sample id. */
	uint GetOutgroupSampleId() const { return m_sampleIO.GetOutgroupSampleId(); }

	/** Get name of outgroup sequences. */
	std::set<std::string> GetOutgroupSeqs() { return m_sampleIO.GetOutgroupSeqs(); };

private:
	/** Read sample data. */
	SampleIO m_sampleIO;

	/** Splits within split system. */
	std::vector<Split> m_splits;
};

#endif

