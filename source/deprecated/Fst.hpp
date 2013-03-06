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

#ifndef _F_ST_
#define _F_ST_

#include "Precompiled.hpp"

#include "SplitSystem.hpp"

/**
 * @class Fst
 * @brief Calculate Fst over a split system.
 */
class Fst
{
public:
	Fst() {}

	void Run(SplitSystem& splitSystem, Matrix& dissMatrix);

private:
	uint PairwiseDifference(const std::string& seq1, const std::string& seq2);
};

#endif