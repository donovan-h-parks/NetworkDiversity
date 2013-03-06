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

#ifndef _SAMPLE_IO_
#define _SAMPLE_IO_

#include "Precompiled.hpp"

/**
 * @brief Read file indicating number of times each sequence is found in a sample.
 */
class SampleIO
{
public:		
	/** Constructor. */
	SampleIO();

	/** Destructor. */
	~SampleIO();

	/**
	* @brief Open sample file.
	*
	* @param filename Path to sample file.
	* @return True if file opened successfully, else false.
	*/
	bool Read(const std::string& filename);

	/** Get number of samples. */
	uint GetNumSamples() const { return m_sampleNames.size(); }

	/** Get sample name. */
	std::string GetSampleName(uint index) const { return m_sampleNames.at(index); }

	/** Get number of sequences (including outgroup and missing sequences). */
	uint GetNumSeqs() const { return m_seqNameToId.size() + m_removedSeqIds.size(); }

	/** Get number of ingroup sequences (excludes outgroup and missing sequences). */
	uint GetNumIngroupSeqs() const { return m_seqNameToId.size(); }

	/** Get sequence id. */
	bool GetSeqId(const std::string& name, uint& seqId);

	/** Check for sequence with the specified name. */
	bool IsSeq(const std::string& name);

	/** Get count data for specified sample. */
	void GetData(uint index, std::vector<double>& count, double& totalNumSeq);

	/** Get name of outgroup sequences. */
	std::set<std::string> GetOutgroupSeqs() { return m_outgroupSeqs; };

	/** Check if there is an outgroup. */
	bool IsOutgroup() const { return m_bOutgroup; }

	/** Get outgroup sample id. */
	uint GetOutgroupSampleId() const { return m_outgroupIndex; }

	/** Check for sequences in sample file not in the phylogeny. */
	void CheckForMissingSeqs(const std::set<std::string>& seqsInPhylogeny);

private:
	/** Determine which, if any, sequences belong to the outgroup. */
	void DetermineOutgroupSeqs();

	/** Remove sequences with the specified ids. */
	void RemoveSeqs(const std::set<uint>& seqIdsToRemove);

private:
	/** File stream. */
	std::ifstream m_file;

	/** Start of each sample in sample count file. */
	std::vector<std::streampos> m_sampleStreamPos;

	/** Name of sequences (only for ingroup sequences). */
	std::map<std::string, uint> m_seqNameToId;

	/** Name of samples. */
	std::vector<std::string> m_sampleNames;

	/** Name of outgroup sequences. */
	std::set<std::string> m_outgroupSeqs;

	/** Seq ids of outgroup sequences. */
	std::set<uint> m_outgroupSeqIds;

	/** Sequences remove from consideration. */
	std::set<uint> m_removedSeqIds;

	/** Temporary buffer for reading sample data. */
	char* m_buffer;

	/** Flag indicating if there is an outgroup sample. Must be labelled 'outgroup' or 'Outgroup'. */
	bool m_bOutgroup;

	/** Index of outgroup sample. */
	uint m_outgroupIndex;
};

#endif

