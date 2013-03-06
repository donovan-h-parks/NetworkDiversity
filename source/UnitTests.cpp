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

#include "UnitTests.hpp"

#include "SplitSystem.hpp"
#include "DiversityCalculator.hpp"

std::string gTempDissFile = "../unit-tests/unit-test.tmp.diss";

bool UnitTests::Execute()
{
	std::cout << "Running unit tests:" << std::endl;

	std::cout << "  Simple tree (implicitly rooted tree, qualitative measures)... ";
	if(!SimpleTreeQual("", "../unit-tests/SimpleTree_ImplicitlyRooted.tre", "../unit-tests/SimpleTree_ImplicitlyRooted.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Simple tree (explicitly rooted tree, qualitative measures)... ";
	if(!SimpleTreeQual("", "../unit-tests/SimpleTree_ExplicitlyRooted.tre", "../unit-tests/SimpleTree_ExplicitlyRooted.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Simple tree (explicitly rooted tree, quantitative measures)... ";
	if(!SimpleTreeQuan("", "../unit-tests/SimpleTree_ExplicitlyRooted.tre", "../unit-tests/SimpleTree_ExplicitlyRooted.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Simple Tree (unrooted split system, quantitative measures)... ";
	if(!SimpleTreeUnrooted("../unit-tests/SimpleTree_Unrooted.nex", "", "../unit-tests/SimpleTree_Unrooted.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Simple Tree (rooted split system, qualitative measures)... ";
	if(!SimpleTreeQual("../unit-tests/SimpleTree_Rooted.nex", "", "../unit-tests/SimpleTree_Rooted.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in tree file... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeQual("", "../unit-tests/SimpleTree_ImplicitlyRooted.tre", "../unit-tests/SimpleTree_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in unrooted split system... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeUnrooted("../unit-tests/SimpleTree_Unrooted.nex", "", "../unit-tests/SimpleTree_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in rooted split system... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeQual("../unit-tests/SimpleTree_Rooted.nex", "", "../unit-tests/SimpleTree_Rooted_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in sample file with rooted tree... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeQual("", "../unit-tests/SimpleTree_ExplicitlyRooted_ExtraSeqs.tre", "../unit-tests/SimpleTree_ExplicitlyRooted_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in sample file and in rooted tree... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeQual("", "../unit-tests/SimpleTree_Rooted_ExtraSeqs.tre", "../unit-tests/SimpleTree_Rooted_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing missing taxa in sample file and in rooted split system... " << std::endl;
	std::cout << "    ";
	if(!SimpleTreeQual("../unit-tests/SimpleTree_Rooted_ExtraSeqs.nex", "", "../unit-tests/SimpleTree_Rooted_ExtraSeqs.env"))
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing multifurcating tree... ";
	if(!Multifurcating())
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	std::cout << "  Testing samples with overlapping sets of sequences... ";
	if(!SharedSeqs())
	{
		std::cout << "failed." << std::endl;
		return false;
	}
	std::cout << "passed." << std::endl;

	return true;
}

bool UnitTests::SimpleTreeQual(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile)
{
	std::vector< std::vector<double> > dissMatrix;

	SplitSystem splitSystem;
	if(!splitSystem.LoadData(nexusFile, newickFile, sampleFile))
		return false;

	// unweighted Bray-Curtis (Sorensen)
	DiversityCalculator uBC(splitSystem, "Bray-Curtis", false);
	if(!uBC.IsGood())
		return false;

	uBC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Canberra
	DiversityCalculator uCanberra(splitSystem, "Canberra", false);
	if(!uCanberra.IsGood())
		return false;

	uCanberra.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted coefficient of similarity
	DiversityCalculator uCS(splitSystem, "CoefficientOfSimilarity", false);
	if(!uCS.IsGood())
		return false;

	uCS.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted CT
	DiversityCalculator uCT(splitSystem, "Complete tree", false);
	if(!uCT.IsGood())
		return false;

	uCT.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Euclidean
	DiversityCalculator uEuclidean(splitSystem, "Euclidean", false);
	if(!uEuclidean.IsGood())
		return false;

	uEuclidean.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][1], sqrt(2.0)))
		return false;

	// unweighted Gower
	DiversityCalculator uGower(splitSystem, "Gower", false);
	if(!uGower.IsGood())
		return false;

	uGower.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted Kulczynski
	DiversityCalculator uKulczynski(splitSystem, "Kulczynski", false);
	if(!uKulczynski.IsGood())
		return false;

	uKulczynski.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 0.5))
		return false;

		// unweighted Kulczynski
	DiversityCalculator uLennonCD(splitSystem, "LCD", false);
	if(!uLennonCD.IsGood())
		return false;

	uLennonCD.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 0.5))
		return false;

	// unweighted Manhattan
	DiversityCalculator uManhattan(splitSystem, "Manhattan", false);
	if(!uManhattan.IsGood())
		return false;

	uManhattan.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted Morisita-Horn
	DiversityCalculator uMH(splitSystem, "Morisita-Horn", false);
	if(!uMH.IsGood())
		return false;

	uMH.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Soergel
	DiversityCalculator uSoergel(splitSystem, "Soergel", false);
	if(!uSoergel.IsGood())
		return false;

	uSoergel.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/3.0))
		return false;

	// unweighted Tamas coefficent
	DiversityCalculator uTC(splitSystem, "Tamas coefficient", false);
	if(!uTC.IsGood())
		return false;

	uTC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted weighted correlation
	DiversityCalculator uWC(splitSystem, "WeightedCorrelation", false);
	if(!uWC.IsGood())
		return false;

	uWC.Dissimilarity(gTempDissFile);	
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1.57735))
		return false;
	if(!Compare(dissMatrix[2][0], 1.57735))
		return false;
	if(!Compare(dissMatrix[2][1], 1.0))
		return false;

	// unweighted Yue-Clayton
	DiversityCalculator uYC(splitSystem, "Yue-Clayton", false);
	if(!uYC.IsGood())
		return false;

	uYC.Dissimilarity(gTempDissFile);	
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/3.0))
		return false;

	return true;
}

bool UnitTests::SimpleTreeQuan(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile)
{
	std::vector< std::vector<double> > dissMatrix;

	SplitSystem splitSystem;
	if(!splitSystem.LoadData(nexusFile, newickFile, sampleFile))
		return false;

	// unweighted Bray-Curtis (Sorensen)
	DiversityCalculator uBC(splitSystem, "Bray-Curtis", true);
	if(!uBC.IsGood())
		return false;

	uBC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Canberra
	DiversityCalculator uCanberra(splitSystem, "Canberra", true);
	if(!uCanberra.IsGood())
		return false;

	uCanberra.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted coefficient of similarity
	DiversityCalculator uCS(splitSystem, "CoefficientOfSimilarity", true);
	if(!uCS.IsGood())
		return false;

	uCS.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted CT
	DiversityCalculator uCT(splitSystem, "Complete tree", true);
	if(!uCT.IsGood())
		return false;

	uCT.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Euclidean
	DiversityCalculator uEuclidean(splitSystem, "Euclidean", true);
	if(!uEuclidean.IsGood())
		return false;

	uEuclidean.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][1], sqrt(2.0)))
		return false;

	// unweighted Gower
	DiversityCalculator uGower(splitSystem, "Gower", true);
	if(!uGower.IsGood())
		return false;

	uGower.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted Kulczynski
	DiversityCalculator uKulczynski(splitSystem, "Kulczynski", true);
	if(!uKulczynski.IsGood())
		return false;

	uKulczynski.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 0.5))
		return false;

		// unweighted Kulczynski
	DiversityCalculator uLennonCD(splitSystem, "LCD", true);
	if(!uLennonCD.IsGood())
		return false;

	uLennonCD.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 0.5))
		return false;

	// unweighted Manhattan
	DiversityCalculator uManhattan(splitSystem, "Manhattan", true);
	if(!uManhattan.IsGood())
		return false;

	uManhattan.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// unweighted Morisita-Horn
	DiversityCalculator uMH(splitSystem, "Morisita-Horn", true);
	if(!uMH.IsGood())
		return false;

	uMH.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted Soergel
	DiversityCalculator uSoergel(splitSystem, "Soergel", true);
	if(!uSoergel.IsGood())
		return false;

	uSoergel.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/3.0))
		return false;

	// unweighted Tamas coefficent
	DiversityCalculator uTC(splitSystem, "Tamas coefficient", true);
	if(!uTC.IsGood())
		return false;

	uTC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// unweighted weighted correlation
	DiversityCalculator uWC(splitSystem, "WeightedCorrelation", true);
	if(!uWC.IsGood())
		return false;

	uWC.Dissimilarity(gTempDissFile);	
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1.57735))
		return false;
	if(!Compare(dissMatrix[2][0], 1.57735))
		return false;
	if(!Compare(dissMatrix[2][1], 1.0))
		return false;

	// unweighted Yue-Clayton
	DiversityCalculator uYC(splitSystem, "Yue-Clayton", true);
	if(!uYC.IsGood())
		return false;

	uYC.Dissimilarity(gTempDissFile);	
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 1))
		return false;
	if(!Compare(dissMatrix[2][0], 1))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/3.0))
		return false;

	return true;

	return true;
}

bool UnitTests::SimpleTreeUnrooted(const std::string& nexusFile, const std::string& newickFile, const std::string& sampleFile)
{
	std::vector< std::vector<double> > dissMatrix;

	SplitSystem splitSystem;
	if(!splitSystem.LoadData(nexusFile, newickFile, sampleFile))
		return false;

	// weighted CT
	DiversityCalculator CT(splitSystem, "Complete tree", true);
	if(!CT.IsGood())
		return false;

	CT.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][0], 3.0/4.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/4.0))
		return false;

	// weighted Euclidean
	DiversityCalculator Euclidean(splitSystem, "Euclidean", true);
	if(!Euclidean.IsGood())
		return false;

	Euclidean.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][0], sqrt(3.0)))
		return false;
	if(!Compare(dissMatrix[2][1], sqrt(2.0)))
		return false;

	// weighted Gower
	DiversityCalculator Gower(splitSystem, "Gower", true);
	if(!Gower.IsGood())
		return false;

	Gower.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	// weighted Manhattan
	DiversityCalculator Manhattan(splitSystem, "Manhattan", true);
	if(!Manhattan.IsGood())
		return false;

	Manhattan.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 3))
		return false;
	if(!Compare(dissMatrix[2][0], 3))
		return false;
	if(!Compare(dissMatrix[2][1], 2))
		return false;

	return true;
}

bool UnitTests::Multifurcating()
{
	std::vector< std::vector<double> > dissMatrix;

	SplitSystem splitSystem;
	if(!splitSystem.LoadData("", "../unit-tests/Multifurcating.tre", "../unit-tests/Multifurcating.env"))
		return false;

	// unweighted Soergel
	DiversityCalculator uSoergel(splitSystem, "Soergel", false);
	if(!uSoergel.IsGood())
		return false;

	uSoergel.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 0.855070))
		return false;

	// weighted Bray-Curtis
	DiversityCalculator BC(splitSystem, "Bray-Curtis", true);
	if(!BC.IsGood())
		return false;

	BC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 0.75457500))
		return false;

	return true;
}

bool UnitTests::SharedSeqs()
{
	std::vector< std::vector<double> > dissMatrix;

	SplitSystem splitSystem;
	if(!splitSystem.LoadData("", "../unit-tests/SharedSeqs.tre", "../unit-tests/SharedSeqs.env"))
		return false;

	// unweighted Tamas coefficient
	DiversityCalculator uTC(splitSystem, "Tamas coefficient", false);
	if(!uTC.IsGood())
		return false;

	uTC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 8.0/13.0))
		return false;
	if(!Compare(dissMatrix[2][0], 8.0/13.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/13.0))
		return false;

	// unweighted Soergel
	DiversityCalculator uSoergel(splitSystem, "Soergel", false);
	if(!uSoergel.IsGood())
		return false;

	uSoergel.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 8.0/13.0))
		return false;
	if(!Compare(dissMatrix[2][0], 8.0/13.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0/9.0))
		return false;

	// unweighted Canberra
	DiversityCalculator uCanberra(splitSystem, "Canberra", false);
	if(!uCanberra.IsGood())
		return false;

	uCanberra.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 8.0))
		return false;
	if(!Compare(dissMatrix[2][0], 8.0))
		return false;
	if(!Compare(dissMatrix[2][1], 2.0))
		return false;

	// weighted Soergel
	DiversityCalculator Soergel(splitSystem, "Soergel", true);
	if(!Soergel.IsGood())
		return false;

	Soergel.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 0.55))
		return false;
	if(!Compare(dissMatrix[2][0], 0.615385))
		return false;
	if(!Compare(dissMatrix[2][1], 0.16981))
		return false;

	// weighted Bray-Curtis
	DiversityCalculator BC(splitSystem, "Bray-Curtis", true);
	if(!BC.IsGood())
		return false;

	BC.Dissimilarity(gTempDissFile);
	ReadDissMatrix(gTempDissFile, dissMatrix);

	if(!Compare(dissMatrix[1][0], 0.37931000))
		return false;
	if(!Compare(dissMatrix[2][0], 0.444440))
		return false;
	if(!Compare(dissMatrix[2][1], 0.092783500))
		return false;

	return true;
}

bool UnitTests::ReadDissMatrix(const std::string& dissMatrixFile, std::vector< std::vector<double> >& dissMatrix)
{
	dissMatrix.clear();

	// open dissimilarity file
	std::ifstream dissIn(dissMatrixFile.c_str());
	if(!dissIn.is_open())
	{
		std::cerr << "Unable to open dissimilarity matrix file: " << dissMatrixFile << std::endl;
		return false;
	}

	uint size;
	dissIn >> size;

	for(uint i = 0; i < size; ++i)
	{
		std::string name;
		dissIn >> name;

		std::vector<double> row;
		for(uint j = 0; j < i; ++j)
		{
			double v;
			dissIn >> v;
			row.push_back(v);
		}

		dissMatrix.push_back(row);
	}

	return true;
}

bool UnitTests::Compare(double actual, double expected)
{
	return fabs(actual - expected) < 0.00001;
}