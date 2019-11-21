/***************************************************************************
 *   Copyright (C) 2006 by Benjamin Curley                                 *
 *   curley@tc.bham.ac.uk						   					       *
 *   rhodan@blueyonder.co.uk		                                       *
 * 									   									   *
 *   GeometryMaker: Creates simple 12 vectrex geometries for metallic      *
 *		    clusters						                               *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include "CGeometry.h"
#include "CIco.h"
#include "Utlities.h"
#include "cXYZreadwrite.h"

using namespace std;

void IncorrectInput()
{
	cout << "===================================================================" << endl;
	cout << "Usage: ./structure_gen <Cub|Ico|Dec|PentBi|Oct> <shellsize> <element> <outputfile> *optional* <bond length>" << endl;
	cout << "===================================================================" << endl;
	cout << "Example ./structure_gen Cub 4 Au Cub309" << endl;
	cout << endl;
	cout << "This will generate a 4(309 atoms) CubeOctahedron Gold Cluster to file Cub309.xyz" << endl;
	cout << "Be aware that they are not minimised and therefore have no energy" << endl;
	cout << endl;
	cout << "===================================================================" << endl;
	cout << "Usage: ./structure_gen BM|RND <Cub|Ico|Dec|PentBi|Oct> <coreSize> <coreElement> <shellsize> <ShellElement> <outputfile> *optional* <bond length>" << endl;
	cout << "===================================================================" << endl;
	cout << "Example:\n";
	cout << "./structure_gen BM Cub 3 Au 5 Ag Au55Ag506.xyz" << endl;
}

int clearUp(const string filename)
// This will reread all the information and save it with the number of atoms at the top
// Bit of a hack but not a problem as quick quick program
{
        vector <string> str_Vector;
        ifstream inData;
        ofstream outData;

        inData.open( filename.c_str() );
        if ( !inData ) return EXIT_FAILURE;

        int size = 0;
        string line;
        while ( getline(inData, line) )
        {
                size++;
                str_Vector.push_back(line);
        }
        inData.close();

        outData.open(filename.c_str());
        if(!outData) {
                cout << "Cannot open output file.\n";
                return EXIT_FAILURE;
        }

        outData << size-1 << endl;

        for(int i=0;i < str_Vector.size(); i++)
                outData << str_Vector.at(i).c_str() << endl;

        outData.close();
}

int main(int argc, char *argv[])
{
	int shell_num, coreshellnum;
	string shell;
	string coreshell;
	string filename;
	string element;
	string coreelement;
	string geometry;
	string bimettalic;
	fstream output;
	float radius;
	CGeometry *ptr;
	cXYZReadWrite reader;

	switch(argc)
	{
		case 1: //single arguement
			IncorrectInput();
			return EXIT_SUCCESS;
		case 2:
			IncorrectInput();
			return EXIT_SUCCESS;
		case 3:
			IncorrectInput();
			return EXIT_SUCCESS;
		case 4:
			IncorrectInput();
			return EXIT_SUCCESS;
		case 5:
			geometry = argv[1];
			shell    = argv[2];
			element  = argv[3];
			filename = argv[4];
			radius = 2.7;
			break;
		case 6:
			geometry = argv[1];
			shell    = argv[2];
			element  = argv[3];
			filename = argv[4];
			from_string<float>(radius,argv[5]);
			break;
		case 8:
			bimettalic	= argv[1];
			geometry	= argv[2];
			coreshell	= argv[3];
			coreelement	= argv[4];
			shell		= argv[5];
			element 	= argv[6];
			filename	= argv[7];
			radius		= 2.7;
			break;

		case 9:
			bimettalic	= argv[1];
			geometry	= argv[2];
			coreshell	= argv[3];
			coreelement	= argv[4];
			shell		= argv[5];
			element 	= argv[6];
			filename	= argv[7];
			from_string<float>(radius, argv[8]);
			break;
	}

	output.open(filename.c_str(), ios::out);

	if(!output.is_open())
	{
		cout << "Error Writing to file" << endl;
		exit(1);
	}

	if(geometry == "Cub")
		ptr = new CCubGen;
	else if(geometry == "Ico")
		ptr = new CIcoGen;
	else if(geometry == "Dec")
		ptr = new CDecGen;
	else if(geometry == "PentBi")
		ptr = new CPentBi;
	else if (geometry == "Oct")
		ptr = new COctGen;
	else
	{
		cout << "Invalid Geometry exiting please check input" << endl;
		exit(1);
	}
	from_string<int>(shell_num, shell);
	from_string<int>(coreshellnum, coreshell);
	ptr->setRadius(radius); //must occur before generateCoordinates
        ptr->setShells(shell_num); // For TO allows to check if we have reached inversion point

        // This initiates writing to file
	reader.WriteFirst(&output);

	if(bimettalic == "BM")
	{
		//ok we need to first generate the core
		for(int i = 0; i <= coreshellnum; i++)
		{
			ptr->generateSingleShell(i, coreelement);

			//reader.SetClusterInfo(ptr->GetShell());
                        //reader.WriteFile(&output);
			//ptr->clearShell();
		}
		//now we want to make the outer shell starting form coreshell+1 finishing on shell
		for(int i = coreshellnum+1; i <= shell_num; i++)
		{
			ptr->generateSingleShell(i, element);

			//reader.SetClusterInfo(ptr->GetShell());
                        //reader.WriteFile(&output);
			//ptr->clearShell();
		}
	}
	else
	{
        //reader.WriteFirst(&output);
        
       		for(int i = 0; i <= shell_num; i++)
		{
			ptr->generateSingleShell(i, element);
			
			//reader.SetClusterInfo(ptr->GetShell());
			//reader.WriteFile(&output);
			//ptr->clearShell();
		}
		
        //ptr->generateCoordinates(shell_num, "Au");
	}
        // Get cluster information
	reader.SetClusterInfo(ptr->GetShell());
	if(bimettalic == "RND")
	{
		reader.RandomiseClusterElements(coreelement,element,3);
	}
        // Write to output
	reader.WriteFile(&output);
	output.close();
	delete ptr;	
	
	// Time for a Logsdail clear up hack...
	clearUp(filename);
	
	return EXIT_SUCCESS;
}
