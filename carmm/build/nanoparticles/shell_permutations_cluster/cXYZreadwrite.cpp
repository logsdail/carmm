#include "cXYZreadwrite.h"
#include "cError.h"
#include "Utlities.h"

#include <fstream>
#include <iostream>


cXYZReadWrite::cXYZReadWrite()
{
}

cXYZReadWrite::~cXYZReadWrite()
{

}

void cXYZReadWrite::ReadFile(std::fstream* pFile)
{
	m_szLine.clear();
	m_szFileContents.clear();
	//check to see that we have a valid open stream


	//deref pointer and grab total files
	while(!pFile->eof())
	{
		getline(*pFile,m_szLine); //shouldnt expect a 1000+ character on a single line
		m_szFileContents.push_back(m_szLine);
	}
	//done with file might aswell close it.
	pFile->close();

	//for some reason we grab one line to many containing nothign so pop it off
	m_szFileContents.pop_back();

	//now we have all the lines but we best check the number of lines which is atoms+2

	if(!from_string<int>(m_ClusterInfo.size, m_szFileContents[0]))
	{
		throw cError("First Line of inputfile is not an integer");
	}

	if(m_ClusterInfo.size + 2 != m_szFileContents.size())
	{
		throw cError("Number of atoms does not match number found in file\n");
	}

	//If we are still here we can start to parse the file
	m_ClusterInfo.comment = m_szFileContents[1];

	for(size_t i = 2; i < m_szFileContents.size(); i++)
	{
		std::vector<std::string> vszWords;
		atomData tempAtom;
		while(m_szFileContents[i].size() > 0)
		{
			//first find any character that matches " \t\n"
			size_t pos = m_szFileContents[i].find_first_of(" \t");
			if(pos != std::string::npos)
			{
				//Ok found a blank space keep removing blankspaces from the string
				//lets find a nonblank space from this position
				pos = m_szFileContents[i].find_first_not_of(" \t",pos);
				vszWords.push_back(m_szFileContents[i].substr(0,pos));
				m_szFileContents[i] = m_szFileContents[i].substr(pos);
			}
			else // no more blank spaces
			{
				vszWords.push_back(m_szFileContents[i]);
				break;
			}

		}
		//best strip any trailling whitespace from atomtype
		size_t pos = vszWords[0].find_first_of(" \t");
		tempAtom.atomType = vszWords[0].substr(0,pos);
		from_string<double>(tempAtom.x,vszWords[1]);
		from_string<double>(tempAtom.y,vszWords[2]);
		from_string<double>(tempAtom.z,vszWords[3]);
		m_ClusterInfo.atom.push_back(tempAtom);
	}	
}	

void cXYZReadWrite::WriteFile(std::fstream *pFile)
{
	if(pFile->is_open())
	{
		*pFile << m_ClusterInfo;
	}
	else
	{
		throw cError("File stream to write to is not present");
	}
}

//should I make my own copy and point to that instead?
//I must remember that if the cluster structure goes out of
//scope then the pointer becomes invalid.
void cXYZReadWrite::SetClusterInfo(clusterData input)
{
	m_ClusterInfo = input;
}	

//think its better I return a copy rather than a pointer
clusterData cXYZReadWrite::GetClusterInfo()
{
	return m_ClusterInfo;
}
