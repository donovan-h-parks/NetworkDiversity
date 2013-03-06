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

#include "Split.hpp"

Split::Split(uint id, double weight, const std::vector<bool>& split, bool bOutgroupSeqOnLeft, bool bOutgroupSeqOnRight, uint numSeqOnLeft, uint numSeqOnRight)
	: m_rightSeqCount(numSeqOnRight), m_leftSeqCount(numSeqOnLeft), m_bOutgroupSeqOnRight(bOutgroupSeqOnRight),
		m_bOutgroupSeqOnLeft(bOutgroupSeqOnLeft), m_split(split), m_weight(weight), m_id(id)
{

}

std::vector<uint> Split::GetLeftSequenceIds() const
{
	std::vector<uint> seqs(m_leftSeqCount);

	uint index = 0;
	for(uint i = 0; i < m_split.size(); ++i)
	{
		if(m_split[i])
		{
			seqs[index] = i;
			index++;
		}
	}

	return seqs;
}

std::vector<uint> Split::GetRightSequenceIds() const
{
	std::vector<uint> seqs(m_rightSeqCount);

	uint index = 0; 
	for(uint i = 0; i < m_split.size(); ++i)
	{
		if(!m_split[i])
		{
			seqs[index] = i;
			index++;
		}
	}

	return seqs;
}

std::vector<uint> Split::GetSequencesIdsInSmallestBipartition() const
{
	if(m_leftSeqCount < m_rightSeqCount)
		return GetLeftSequenceIds();

	return GetRightSequenceIds();
}

bool Split::IsTrivial() const
{
	// check if left or right side of split contains a single sequence
	return ((m_leftSeqCount == 1 && m_bOutgroupSeqOnRight) || (m_rightSeqCount == 1 && m_bOutgroupSeqOnLeft));
}

int Split::GetSize() const
{
	return std::min<uint>(m_leftSeqCount, m_rightSeqCount);
}
