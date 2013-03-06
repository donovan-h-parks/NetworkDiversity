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

#include "SampleIO.hpp"
#include "Utils.hpp"

std::string TrimStr(const std::string& Src, const std::string& c = " \r\n")
{
	int p2 = Src.find_last_not_of(c);

	if (p2 == std::string::npos) 
		return std::string();

	int p1 = Src.find_first_not_of(c);
	if (p1 == std::string::npos) 
		p1 = 0;

	return Src.substr(p1, (p2-p1)+1);
}

SampleIO::SampleIO(): m_bOutgroup(false), m_outgroupIndex(std::numeric_limits<uint>::max())
{
	m_buffer = NULL;
}

SampleIO::~SampleIO() 
{ 
	if(m_buffer)
		delete[] m_buffer;

	if(m_file.is_open())
		m_file.close(); 
}

bool SampleIO::Read(const std::string& filename)
{
	m_file.open(filename.c_str());
	if(!m_file.is_open())
	{
		std::cerr << "Unable to open sample file: " << filename << std::endl;
		return false;
	}

	// check if file ends with a end-of-line character(s)
	char c;
	bool bEndOfLineTerminator = true;
	m_file.seekg(-1, std::ios::end);
	m_file.read(&c, 1);
	if(c != '\n')
		bEndOfLineTerminator = false;
	m_file.seekg(0, std::ios::beg);

	// parse header line to get order of sequences
	std::string line, token;
	std::getline(m_file, line);
	std::stringstream ss(line);
	uint seqId = 0;
	while(std::getline(ss, token, '\t'))
	{
		if(!token.empty())
		{
			m_seqNameToId[TrimStr(token)] = seqId;
			seqId++;
		}
	}

	// get number of samples and starting index of each sample line
	std::streamsize longestLine = 0;
	do
	{
		if(!line.empty())
		{
			if(!m_sampleStreamPos.empty())
			{
				uint pos = line.find('\t');
				std::string sampleName = line.substr(0, pos);
				if(sampleName == "outgroup" || sampleName == "Outgroup")
				{
					m_bOutgroup = true;
					m_outgroupIndex = m_sampleNames.size();
				}
				else
					m_sampleNames.push_back(sampleName);

				std::streamsize lineLen = m_file.tellg() - m_sampleStreamPos[m_sampleStreamPos.size()-1] + 2; // +2 is for end-of-line character
				if(lineLen > longestLine)
					longestLine = lineLen;
			}

			m_sampleStreamPos.push_back(m_file.tellg());
		}
	} while(std::getline(m_file, line));
	m_file.clear();

	// allocate temporary buffer for reading lines
	m_buffer = new char[longestLine];

	if(!bEndOfLineTerminator)
	{
		int endOfLineLen = 1;
		#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
			endOfLineLen = 2;
		#endif
	
		m_file.seekg(endOfLineLen, std::ios::end);
		m_sampleStreamPos[m_sampleStreamPos.size()-1] = m_file.tellg();
	}

	DetermineOutgroupSeqs();

	return true;
}

void SampleIO::GetData(uint index, std::vector<double>& count, double& totalNumSeq)
{
	// read the ith sample from file
	m_file.clear();
	m_file.seekg(m_sampleStreamPos[index]);
	
	std::streamsize charsInLine = m_sampleStreamPos[index+1] - m_sampleStreamPos[index] - 1;
	m_file.read(m_buffer, charsInLine);
	m_buffer[charsInLine] = 0;
	
	// read sample name
	char* curPos = (char *)memchr(m_buffer, '\t', (size_t)charsInLine);
	++curPos;
	charsInLine -= (curPos - m_buffer);

	// read count data
	totalNumSeq = 0;
	count.clear();
	count.reserve(GetNumIngroupSeqs());
	uint seqId = 0;
	do
	{
		char* tabPos = (char *)memchr(curPos, '\t', charsInLine);
		if(tabPos != NULL)
			tabPos[0] = 0;

		double numSeq = fast_atof(curPos);
		if(m_removedSeqIds.count(seqId) == 0)
		{
			count.push_back(numSeq);
			totalNumSeq += numSeq;		
		}

		charsInLine -= (tabPos - curPos) + 1;
		curPos = tabPos + 1;
		seqId++;
	}while(count.size() != GetNumIngroupSeqs());
}

void SampleIO::DetermineOutgroupSeqs()
{
	if(m_bOutgroup)
	{
		// determine outgroup sequences and remove from ingroup 
		std::vector<double> count;
		double totalNumSeq;
		GetData(m_outgroupIndex, count, totalNumSeq);
		std::map<std::string, uint> ingroupSeqNameToId;
		for(uint seqId = 0; seqId < count.size(); ++seqId)
		{
			if(count[seqId] > 0)
			{
				m_outgroupSeqIds.insert(seqId);

				// find sequence with the specified id
				std::map<std::string, uint>::iterator iter;
				for (iter = m_seqNameToId.begin(); iter != m_seqNameToId.end(); ++iter)
				{
					if(iter->second == seqId)
					{
						m_outgroupSeqs.insert(iter->first);
						break;
					}
				}
			}
		}
	}
}

void SampleIO::CheckForMissingSeqs(const std::set<std::string>& seqsInPhylogeny)
{
	// find sequences in sample file not contained in the tree
	std::set<std::string> seqsMissingInTree;
	std::map<std::string, uint>::iterator iter;
	for(iter = m_seqNameToId.begin(); iter != m_seqNameToId.end(); ++iter)
	{
		std::string name = iter->first;

		if(seqsInPhylogeny.count(name) == 0)
			seqsMissingInTree.insert(name);
	}

	// report missing sequences
	if(seqsMissingInTree.size() > 0)
	{
		std::cout << "(Warning) The following taxa are in your sample file, but not your phylogeny file:" << std::endl;

		std::set<std::string>::iterator iter;
		for(iter = seqsMissingInTree.begin(); iter != seqsMissingInTree.end(); ++iter)
			std::cout << *iter << std::endl;
	}

	// remove missing sequences from sample file
	std::set<uint> seqIdsToRemove = m_outgroupSeqIds;
	std::set<std::string>::iterator it;
	for(it = seqsMissingInTree.begin(); it != seqsMissingInTree.end(); ++it)
	{
		uint seqId;
		GetSeqId(*it, seqId);
		seqIdsToRemove.insert(seqId);
	}

	RemoveSeqs(seqIdsToRemove);
}

void SampleIO::RemoveSeqs(const std::set<uint>& seqIdsToRemove)
{
	// remove specified sequences and ensure sample data is return
	// with these sequences being ignored
	//
	// Note: this function should only be called once will all 
	// sequnces to be removed.

	std::map<uint, std::string> seqIdToSeqName;
	std::map<std::string, uint>::iterator iter;
	for(iter = m_seqNameToId.begin(); iter != m_seqNameToId.end(); ++iter)
		seqIdToSeqName[iter->second] = iter->first;

	std::map<std::string, uint> ingroupSeqNameToId;
	uint removedSeqCount = 0;
	std::map<uint, std::string>::iterator iter2;
	for(iter2 = seqIdToSeqName.begin(); iter2 != seqIdToSeqName.end(); ++iter2)
	{
		if(seqIdsToRemove.count(iter2->first) == 0)
			ingroupSeqNameToId[iter2->second] = iter2->first - removedSeqCount;
		else
		{
			m_removedSeqIds.insert(iter2->first);
			removedSeqCount++;
		}
	}

	m_seqNameToId = ingroupSeqNameToId;
}

bool SampleIO::GetSeqId(const std::string& name, uint& seqId) 
{ 
	std::map<std::string, uint>::iterator it;
	it = m_seqNameToId.find(name);
	if(it != m_seqNameToId.end())
	{
		seqId = it->second;
		return true;
	}

	seqId = 0;
	return false;
}
