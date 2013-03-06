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

#ifndef _NEXUS_IO_
#define _NEXUS_IO_

#include "Precompiled.hpp"

class SplitSystem;

/**
 * @brief Populate a split system from a Nexus file containing a TAXA and either a TREES or SPLITS blocks.
 */
class NexusIO
{
public:
	/** Constructor. */
	NexusIO() {}

	/** Destructor. */
	~NexusIO() {}

	/**
	* @brief Read a split network from a file.
	*
	* @param splitSystem Split system to populate from file.
	* @param filename The file path.
	* @return True if file parsed successfully, otherwise false.
	*/
	bool Read(SplitSystem *const splitSystem, const std::string& filename);

private:
	/**
	* @brief Read TAXA block.
	*
	* @param splitSystem Split system to populate from file.
	* @param textStream Stream for reading Nexus input file.
	* @return True if block parsed successfully, otherwise false.
	*/
	bool ReadTaxaBlock(SplitSystem *const splitSystem, std::ifstream& textStream);

	/**
	* @brief Read CHARACTERS block.
	*
	* @param splitSystem Split system to populate from file.
	* @param textStream Stream for reading Nexus input file.
	* @return True if block parsed successfully, otherwise false.
	*/
	bool ReadCharactersBlock(SplitSystem *const splitSystem, std::ifstream& textStream);

	/**
	* @brief Read TREES block.
	*
	* @param splitSystem Split system to populate from file.
	* @param textStream Stream for reading Nexus input file.
	* @return True if block parsed successfully, otherwise false.
	*/
	bool ReadTreesBlock(SplitSystem *const splitSystem, std::ifstream& textStream);

	/**
	* @brief Read SPLITS block.
	*
	* @param splitSystem Split system to populate from file.
	* @param textStream Stream for reading Nexus input file.
	* @return True if block parsed successfully, otherwise false.
	*/
	bool ReadSplitsBlock(SplitSystem *const splitSystem, std::ifstream& textStream);

private:
	/** Map nexus id to sequence name. */
	std::map<uint, std::string> m_nexusIdToName;
};


#endif


