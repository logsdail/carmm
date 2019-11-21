/*

	Program Flow:
	Initalise Minimiser
	Input Coordinates.
	Initalise Permutations.
	retrive permutations.
    Alter structure based on permutations.

*/

#include <fstream>
#include <string>
#include "cPermute.h"		  //include permutation class.
#include "cgupta_minimiser.h" //include guptaminimiser
#include "cXYZreadwrite.h" 	  //include xyzreadwrite class
#include "Utlities.h"

using namespace std;

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
        // Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

/*  Inputfile of the form
    potential
    structure(must have dummy atom positions)
    ElementA ElementB
    int int int int int int ... (for each atom that represents start of subshell)
    int int int int int int ... (for each atom that represents end of subshell)
    prefix (what is added at the begining of each filename).
    minimise or not
    if minimise, cutoff or not
*/

int main(int argc, char **argv)
{
	//get input 
	string potential_filename;
	string structure_filename;
	string elementA, elementB;
	string inputbuffer;
	string prefix;

	bool bCutoff = false;
	bool bMinimise = false;

	vector<std::string> tokens;

	//hold first atom pos and last atom pos of a subshell
	vector<int> vStartPositions;
	vector<int> vEndPositions;

	//hold number of atoms in each shell
	vector<int> vAtomsInShell;

	fstream inputStream;
	inputStream.open(argv[1], std::ios::in);
	if(!inputStream.is_open())
	{
		std::cout << "Unable to find input file" << std::endl;
		exit(1);
	}

	getline(inputStream, inputbuffer); //pot
	Tokenize(inputbuffer, tokens);
	getline(inputStream, inputbuffer); //strut
	Tokenize(inputbuffer, tokens);
	getline(inputStream, inputbuffer); //elements
	Tokenize(inputbuffer, tokens);

	potential_filename = tokens[0];
	structure_filename = tokens[1];
	elementA = tokens[2];
	elementB = tokens[3];
	tokens.clear();

	getline(inputStream, inputbuffer);
	Tokenize(inputbuffer, tokens);
	for(int i = 0; i < tokens.size(); i++)
	{
		int temp;
		from_string(temp, tokens[i]);
		temp -= 1; //this as people use 1..n rather than 0..n-1
		vStartPositions.push_back(temp);
	}

	tokens.clear();
	getline(inputStream, inputbuffer);
	Tokenize(inputbuffer, tokens);
	for(int i = 0; i < tokens.size(); i++)
	{
		int temp;
		from_string(temp, tokens[i]);
		temp -= 1;
		vEndPositions.push_back(temp);
	}
	tokens.clear();
	getline(inputStream, inputbuffer); // Prefix
	Tokenize(inputbuffer, tokens);
        getline(inputStream, inputbuffer); // Minimise?
        Tokenize(inputbuffer, tokens);
	getline(inputStream, inputbuffer); // Cutoff?
        Tokenize(inputbuffer, tokens);

	prefix = tokens[0];
	if (tokens[1] == "minimise") 
	{
		bMinimise = true;
		if (tokens[2] == "cutoff") bCutoff = true;
	}
	
	//Work out number of atoms in each shell (Reliant on start/end positions matching)
        for (int i = 0; i < vStartPositions.size(); i++)
        {
                int temp;
		temp = vEndPositions[i] + 1;
                temp -= vStartPositions[i];
                vAtomsInShell.push_back(temp);
        }

	//Ok all stuff loaded, not neatly but done
	//create minimiser
	CGupta_minimiser minimiser(potential_filename,bCutoff);
	//create permutator set it to size of number of subshells
	cPermute Permutations(vStartPositions.size());
	Permutations.BeginPermutations();

	std::fstream structureStreamIn;
	std::fstream structureStreamOut;
	structureStreamIn.open(structure_filename.c_str(),std::ios::in);
	cXYZReadWrite structureIn;
	cXYZReadWrite structureOut;

	structureIn.ReadFile(&structureStreamIn); //read structure file in

	std::list<std::vector<int> >  lvPermutations;
	lvPermutations = Permutations.GetPermuationList();

	std::list<std::vector<int> >::iterator iter;
	//for each permutation
	for(iter = lvPermutations.begin(); iter != lvPermutations.end(); iter++)
	{
		int countElementA = 0;
		int countElementB = 0;
		clusterData temp = structureIn.GetClusterInfo();
		std::string outputfile;
		std::string binaryTemp;
		outputfile = prefix;

		for(int i = 0; i < iter->size(); i++)
		{
			binaryTemp += to_string((*iter)[i]);
			//change structure
			if((*iter)[i] == 0)
			{
				countElementA += vAtomsInShell[i];
				for(int currentPos = vStartPositions[i]; currentPos <= vEndPositions[i]; currentPos++)
				{
					temp.atom[currentPos].atomType = elementA;
				}	
			}
			if((*iter)[i] == 1)
			{
				countElementB += vAtomsInShell[i];
				for(int currentPos = vStartPositions[i]; currentPos <= vEndPositions[i]; currentPos++)
				{
					temp.atom[currentPos].atomType = elementB;
				}	
			}
		}
		//ok temp now holds our new homotop arrangement, need to minimise and save structure
		minimiser.SetClusterData(temp);
		//Editing of filename to include compositive structure as well as binary of shells
		outputfile += elementA + to_string(countElementA);
		outputfile += elementB + to_string(countElementB);
		outputfile += "_" + binaryTemp;
		outputfile += ".xyz";
		structureStreamOut.open(outputfile.c_str(),std::ios::out);
		
                // First Im going to edit it to not minimise structures
                // Then we'll look at adding a boolean operator to choose minimise or not
                // structureOut.SetClusterInfo(minimiser.GetClusterData());
                // structureOut.SetClusterInfo(minimiser.GetMinimisedData()); //this minimises it

                // Sorted - so this will choose to minimise or not              
                if (bMinimise) structureOut.SetClusterInfo(minimiser.GetMinimisedData());
		else structureOut.SetClusterInfo(minimiser.GetClusterData());

		structureOut.WriteFile(&structureStreamOut);
		structureStreamOut.close();
		std::cout << outputfile << "\t" << structureOut.GetClusterInfo().potentialEnergy << std::endl;
	}
}

