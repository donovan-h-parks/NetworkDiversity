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

#ifndef _UNIT_TESTS_

#include "Precompiled.hpp"

/**
 * @brief Execute unit tests.
 */
class UnitTests
{
public:		
	bool Execute();

private:
	bool ReadDissMatrix(const std::string& dissMatrixFile, std::vector< std::vector<double> >& dissMatrix);
	bool Compare(double actual, double expected);

	/** Test several variants of a simple tree. Ground truth determined by hand.
	 *   a) implicitly rooted tree 
	 *   b) explicitly rooted tree
	 *   c) unrooted split system
	 *   d) explicitly rooted split system
	 *   e) missing sequences from sample file and/or phylogeny
	*/
	bool SimpleTreeQual(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile);
	bool SimpleTreeQuan(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile); 	
	bool SimpleTreeUnrooted(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile); 	

	/** Test multifuricating tree. Ground truth determined by Chameleon and Fast UniFrac. */
	bool Multifurcating();

	/** Test tree with shared sequences. Ground truth determined by Chameleon and Fast UniFrac. */
	bool SharedSeqs();
};

#endif

