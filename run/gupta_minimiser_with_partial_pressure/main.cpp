#include "cXYZreadwrite.h"
#include "cgupta_minimiser.h"
#include <iostream>
#include <fstream>
#include <string>

#include "cError.h"
#include "cmorse_minimiser.h"

using namespace std;

int main(int argc, char *argv[])
{	
	bool bMinimise = true;
	bool bCutoff = false;
	bool bVerbose = false;

	if(argc != 2 && argc != 3)
	{
		//we either dont have an inputfile or to many arguments
		std::cout << "USAGE: ./Minimiser <inputfilename> (verbose)" << std::endl;
		std::cout << std::endl;
		std::cout << "Inputfile name must contain in the correct order:" << std::endl;
		std::cout << "potential.in" << std::endl;
		std::cout << "Struture.xyz" << std::endl;
		std::cout << "minimisedStructure.xyz" << std::endl;
		std::cout << "remin | calc" << std::endl;
		std::cout << "(cutoff)" << std::endl;
		std::cout << std::endl;
                std::cout << "Inputs in brackets are optional" << std::endl;
		std::cout << "If used cutoff requires cutoff_parameters file" << std::endl;
		return EXIT_FAILURE;
	}

	if (argc == 3)
	{
		std::string verbose = argv[2];
		if (verbose == "verbose")
		{
			bVerbose = true;
		}
		else
		{
			std::cout << "Unknown command: " << verbose << ". Please consult options." << endl;
			return EXIT_FAILURE;
		}
	}

	std::fstream input;
	std::string inputname = argv[1];
	input.open(inputname.c_str(), std::ios::in);
	if(!input.is_open())
	{
		throw cError ("Unable to find input file");
	}

	std::string potential;
	std::string structure;
	std::string outputstructure;
	std::string option;
	std::string cutoff;

	//get input files contents
	input >> potential;
	input >> structure;
	input >> outputstructure;
	input >> option;
	input >> cutoff;

	if(option == "remin")
		bMinimise = true;
	if(option == "calc")
		bMinimise = false;

        if(cutoff == "cutoff")
                bCutoff = true;

	//create the xyz IO class
	cXYZReadWrite reader;
	//create the minimiser
	CGupta_minimiser theMinimiser(potential, bCutoff, bVerbose);
	std::fstream file_op;
	std::fstream file_ot;

	//open both streams
	file_op.open(structure.c_str(), ios::in);
	if(!file_op.is_open())
	{
		std::string errorText = "Structure File: " + structure + " does not exist!";
		throw cError(errorText.c_str());
	}
	file_ot.open(outputstructure.c_str(), ios::out);
	if(!file_ot.is_open())
	{
		throw cError("Unable to create output file!");
	}

	//read in file
	reader.ReadFile(&file_op);

	//pass the minimiser the structure info
	theMinimiser.SetClusterData(reader.GetClusterInfo());
	if(bMinimise)
	{
		//update the reader with new minimised structure info
		reader.SetClusterInfo(theMinimiser.GetMinimisedData());
	}
	else
	{
		clusterData temp;
		temp = reader.GetClusterInfo();
		temp.potentialEnergy = theMinimiser.GetEnergy();
		reader.SetClusterInfo(temp);
	}
	//write the file
	reader.WriteFile(&file_ot);

	return EXIT_SUCCESS;
}
