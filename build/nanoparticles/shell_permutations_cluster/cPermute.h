/*
	Simple class to produce all possible permutaitons of a fixed size set.

	Should be easly expanded to permute more than 2 choice at each position

*/


#include <vector>
#include <iostream>
#include <list>

using std::vector; using std::list;
using std::cout;   using std::endl;

class cPermute
{
public:
	cPermute(int size) : m_iSize(size) {;}
	~cPermute();
	void BeginPermutations();
	void WritePremutations();
	list<vector<int> > GetPermuationList(){return m_lvPermutations;}
private:
	void permute(int index);
	int m_iSize;
	vector<int> m_vPermutationArray;
	list<vector<int> > m_lvPermutations;
};
