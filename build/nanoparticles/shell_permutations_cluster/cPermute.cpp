#include "cPermute.h"

cPermute::~cPermute()
{


}

void cPermute::BeginPermutations()
{
	//make sure we starting a fresh
	m_vPermutationArray.clear();
	m_lvPermutations.clear();
	m_vPermutationArray.resize(m_iSize);
	permute(m_iSize-1);
}


//Debug Method check that permutations are stored.
void cPermute::WritePremutations()
{
	list<vector<int> >::iterator iter;
	for(iter = m_lvPermutations.begin(); iter != m_lvPermutations.end(); iter++)
	{
		int vecSize = iter->size();
		for(int i = 0; i < vecSize; i++)
		{
			cout << (*iter)[i];
		}
		cout << endl;
	}
}



//TODO: This method can be refactors and made more generic
//Avoidable code repetition.
void cPermute::permute(int index)
{

	if( index == 0 )
	{
		//create temp vector
		vector<int> temp;

		m_vPermutationArray[index] = 0;
		for(int i = 0; i < m_iSize; i++ )
		{
			temp.push_back(m_vPermutationArray[i]);
		}
		m_lvPermutations.push_back(temp);


		//clear our temp vector
		temp.clear();
		m_vPermutationArray[index] = 1;
		for(int i = 0; i < m_iSize; i++ )
		{
			temp.push_back(m_vPermutationArray[i]);
		}
		m_lvPermutations.push_back(temp);
	}
	else
	{
		m_vPermutationArray[index] = 0;
		permute((index -1));

		m_vPermutationArray[index] = 1;
		permute((index -1));
	}

}
