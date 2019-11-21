#include "Struct.h"

ostream& operator<<(ostream& os, const clusterData& c)
{
	//os << "  " << c.size << std::endl;
	//os << std::setiosflags (std::ios::left|std::ios::fixed) << std::setprecision(8) << "   " << c.potentialEnergy << std::endl;

	for(int i = 0; i < c.size; i++)
	{
		os << setiosflags( std::ios::right|std::ios::fixed ) << std::setprecision(12) << c.atom[i].atomType << "  " 
		<< std::setw(16) << c.atom[i].x << "  "
		<< std::setw(16) << c.atom[i].y << "  " 
		<< std::setw(16) << c.atom[i].z ;
		os << std::endl;
	}
	return os; 
}
