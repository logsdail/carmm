#pragma once

#ifndef __CXYZREADER_H__
#define __CXYZREADER_H__

#include <sstream>
#include <iostream>
#include <fstream> //stl file handling
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>

#include "Struct.h"

//Ulitmately this class will use a generic parser which tokenises input
//For now it will be all declared in here

/*General format of XYZ file
	firstline contains the number of atoms in the file
	secondline contains a comment, generally the energy of the structure is placed here
	then each subsequent line has the following form
	String double double double \n
	where string is atom type double are x y z respectively and then the line ends
	Non standard files may contain further columns so after 4 columns non standard data
	may exist.
*/

class cXYZReadWrite
{
public:
	cXYZReadWrite();
	~cXYZReadWrite();

	void ReadFile(std::fstream *pFile); 	 //read the file
	void WriteFirst(std::fstream *pFile);     //write the initial data
	void WriteFile(std::fstream *pFile);     //write the file
	clusterData GetClusterInfo();			 //returns a copy of clusterData
	void SetClusterInfo(clusterData input); //

	//Terrible as this doesn't belong here but hey its a hack
	void RandomiseClusterElements(std::string elementA, std::string elementB, int seed);
private:
	clusterData m_ClusterInfo;
	std::vector<std::string> m_szFileContents; //hold the contents of a file by line
	std::string m_szLine;
};


#endif
