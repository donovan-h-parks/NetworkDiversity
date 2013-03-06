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

#ifndef _NEWICK_IO_
#define _NEWICK_IO_

#include "Precompiled.hpp"

#include "Tree.hpp"
#include "Node.hpp"

class SplitSystem;

/**
 * @brief Populate a tree from a Newick formatted tree file.
 */
class NewickIO
{
public:		
	/** Constructor. */
	NewickIO() {}

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
	* @param Build tree structure from a file containing a Newick formatted tree.
	*
	* @param tree Tree to populate from file.
	* @param filename Path to file containing Newick formatted tree.
	* @return True if tree built successfully, else false if reading file failed.
	*/
	bool Read(Tree<Node>& tree, const std::string& filename);

	/**
	* @brief Parse a string in Newick format and convert it into a tree.
	*
	* @param tree Tree to populate.
	* @param newickStr String containing tree in Newick format.
	* @return True if structures populated successfully, else false if reading file failed.
	*/
	bool ParseNewickString(Tree<Node>& tree, const std::string& newickStr);

	/**
	* @brief Parse Newick string information about a node.
	*
	* @param node Node to associate information with.
	* @param nodeInfo Node data obtained from Newick string.
	*/
	void ParseNodeInfo(Node* node, std::string& nodeInfo);
};

#endif

