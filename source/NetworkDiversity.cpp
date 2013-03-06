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

#include "getopt_pp.hpp"

#include "SplitSystem.hpp"
#include "DiversityCalculator.hpp"

#include "UnitTests.hpp"

bool ParseCommandLine(int argc, char* argv[], std::string& calculator, std::string& nexusFile, 
												std::string& newickFile, std::string& sampleFile, std::string& outputFile, 
												bool& bWeighted, bool& bCount, uint& maxDataVecs, bool& bVerbose)
{
	bool bShowHelp;
	bool bUnitTests;
	bool bShowCalc;
	std::string maxDataVecsStr;
	GetOpt::GetOpt_pp opts(argc, argv);
	opts >> GetOpt::OptionPresent('h', "help", bShowHelp);
	opts >> GetOpt::OptionPresent('l', "list-calc", bShowCalc);
	opts >> GetOpt::OptionPresent('v', "verbose", bVerbose);
	opts >> GetOpt::OptionPresent('u', "unit-tests", bUnitTests);
	opts >> GetOpt::Option('c', "calculator", calculator);
	opts >> GetOpt::Option('n', "nexus-file", nexusFile);
	opts >> GetOpt::Option('t', "newick-file", newickFile);
	opts >> GetOpt::Option('s', "sample-file", sampleFile);
	opts >> GetOpt::Option('o', "output-file", outputFile);
	opts >> GetOpt::Option('x', "max-data-vecs", maxDataVecsStr, "1000");
	opts >> GetOpt::OptionPresent('w', "weighted", bWeighted);
	opts >> GetOpt::OptionPresent('y', "count", bCount);

	maxDataVecs = atoi(maxDataVecsStr.c_str());

	if(bShowHelp || argc <= 1) 
	{		
		std::cout << std::endl;
		std::cout << "Network Diversity v1.0.0 (April 10, 2012)" << std::endl;
		std::cout << "  by Donovan Parks (parks@cs.dal.ca) and Rob Beiko (beiko@cs.dal.ca)" << std::endl;
		std::cout << std::endl;
		std::cout << " Usage: " << opts.app_name() << " -n <nexus file> -s <sample file> -o <output file>" << std::endl;
		std::cout << "        " << opts.app_name() << " -t <tree file> -s <sample file> -o <output file>" << std::endl;
		std::cout << "  -h, --help           Produce help message." << std::endl;
		std::cout << "  -l, --list-calc      List all supported calculators." << std::endl;
		std::cout << "  -u, --unit-tests     Execute unit tests." << std::endl;
		std::cout << std::endl;
		std::cout << "  -c, --calculator     Beta-diversity calculator to use (e.g., Manhattan or uUF)." << std::endl;
		std::cout << "  -w, --weighted       Indicates if sequence abundance data should be used." << std::endl;
		std::cout << "  -y, --count          Use count data as opposed to relative proportions." << std::endl;
		std::cout << std::endl;
		std::cout << "  -n, --nexus-file     Nexus input file (must contain a Taxa and Splits block, root split system with Outgroup taxa)." << std::endl;
		std::cout << "  -t, --newick-file    Newick input file (tree treated as implicitly rooted)." << std::endl;
		std::cout << "  -s, --sample-file    Sample file indicating number of times each sequences is found in a sample." << std::endl;
		std::cout << "  -o, --output-file    Output file." << std::endl;
		std::cout << std::endl;
		std::cout << "  -x, --max-data-vecs  Maximum number of samples to have in memory at once (default = 1000)." << std::endl;
		std::cout << std::endl;
		std::cout << "  -v, --verbose        Provide additional information on program execution." << std::endl;
							
    return false;
  }
	else if(bShowCalc)
	{
		std::cout << std::endl;
		std::cout << "Qualitative measures requiring a rooted split system:" << std::endl;
		std::cout << "  Bray-Curtis (aka: Sorensen, PhyloSor, Dice's index, pairwise Whittaker)" << std::endl;
		std::cout << "  Canberra" << std::endl;
		std::cout << "  Coefficient of similarity (CS)" << std::endl;
		std::cout << "  Euclidean" << std::endl;
		std::cout << "  Gower" << std::endl;
		std::cout << "  Kulczynski (aka: Kulczynski-Cody, Sokal-Sneath)" << std::endl;
		std::cout << "  Lennon compositional difference (LCD)" << std::endl;		
		std::cout << "  Manhattan (aka: Hamming distance)" << std::endl;
		std::cout << "  Soergel (aka: unweighted UniFrac, Jaccard)" << std::endl;
		std::cout << "  Tamas coefficient (TC) (aka: simple matching coefficent)" << std::endl;
		std::cout << "  Weighted correlation (WC)" << std::endl;
		std::cout << std::endl;
		std::cout << "Quantitative measures requiring a rooted split system (use -w flag):" << std::endl;
		std::cout << "  Bray-Curtis (aka: normalized weighted UniFrac, percentage difference)" << std::endl;
		std::cout << "  Canberra" << std::endl;
		std::cout << "  Coefficient of similarity (CS)" << std::endl;
		std::cout << "  Complete tree" << std::endl;
		std::cout << "  Euclidean" << std::endl;
		std::cout << "  Gower" << std::endl;
		std::cout << "  Kulczynski" << std::endl;
		std::cout << "  Lennon compositional difference (LCD)" << std::endl;		
		std::cout << "  Manhattan (aka: weighted UniFrac)" << std::endl;
		std::cout << "  Morisita-Horn" << std::endl;
		std::cout << "  Soergel (aka: Ruzicka, Marczewski-Steinhaus, percentage remoteness)" << std::endl;
		std::cout << "  Tamas coefficient (TC) (aka: simple matching coefficent)" << std::endl;
		std::cout << "  Weighted correlation (WC)" << std::endl;
		std::cout << "  Yue-Clayton (aka: similarity ratio)" << std::endl;
		std::cout << "Quantitative measures requiring a rooted split system (use -w flag):" << std::endl;
		std::cout << "  Complete tree" << std::endl;
		std::cout << "  Euclidean" << std::endl;
		std::cout << "  Gower" << std::endl;
		std::cout << "  Manhattan (aka: weighted UniFrac)" << std::endl;
		std::cout << std::endl;

		return false;
	}
	else if(bUnitTests)
	{
		std::cout << std::endl;

		UnitTests unitTests;
		if(unitTests.Execute())
		{
			std::cout << std::endl;
			std::cout << "Passed all unit tests." << std::endl;
		}
		else
		{
			std::cout << std::endl;
			std::cerr << "Failed unit test." << std::endl;
		}

		return false;
	}

	if(!nexusFile.empty() && !newickFile.empty())
	{
		std::cerr << "Specify either a Nexus (-n) or Newick (-t) file." << std::endl;
		return false;
	}

	if(calculator == "Normalized weighted UniFrac" || calculator == "NWU" || calculator == "NormalizedWeightedUniFrac" || calculator == "Normalized Weighted UniFrac")
	{
		calculator = "Bray-Curtis";
	}

	return true;
}

int main(int argc, char* argv[])
{
	std::clock_t timeStart = std::clock();

	// Parse command line arguments
	std::string calculator;
	std::string nexusFile;
	std::string newickFile;
	std::string sampleFile;
	std::string outputFile;
	bool bWeighted;
	bool bCount;
	bool bVerbose;
	uint maxDataVecs;
	if(!ParseCommandLine(argc, argv, calculator, nexusFile, newickFile, sampleFile, outputFile, bWeighted, bCount, maxDataVecs, bVerbose))
		return 0;

	// create split system
	SplitSystem splitSystem;
	if(!splitSystem.LoadData(nexusFile, newickFile, sampleFile, bVerbose))
	{
		std::cerr << "Failed to load data.";
		return -1;
	}

	// run beta-diversity measures
	DiversityCalculator diversityCalc(splitSystem, calculator, bWeighted, bCount, maxDataVecs, bVerbose);
	if(!diversityCalc.IsGood())
		return -1;

	diversityCalc.Dissimilarity(outputFile);

	std::clock_t timeEnd = std::clock();

	if(bVerbose)
	{
		std::cout << "Total running time: " << ( timeEnd - timeStart ) / (double)CLOCKS_PER_SEC << " s" << std::endl; 
		std::cout << "Done." << std::endl;
	}

	
	return 0;
}

