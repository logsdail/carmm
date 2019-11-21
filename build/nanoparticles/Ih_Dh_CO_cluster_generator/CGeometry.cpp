/***************************************************************************
 *   Copyright (C) 2006 by Benjamin Curley                                 *
 *   curley@tc.bham.ac.uk						  						   *	
 *   rhodan@blueyonder.co.uk		                                       *
 * 									 									   *
 *   GeometryMaker: Creates simple 12 vectrex geometries for metallic      *
 *		    clusters						  							   *
 ***************************************************************************/

#include "CGeometry.h"
#include "Utlities.h"

using namespace std;


//Class Definitions
CGeometry::CGeometry()
{
}

CGeometry::~CGeometry()
{

}

void CGeometry::clearShell()
{
     m_Cluster.atom.clear();
     m_Cluster.size = 0;
     m_Cluster.comment.clear();
}

void CGeometry::add_coordinate(vec3d temp, std::string element)
{
	atomData tempAtom;
	m_Cluster.size++; //every time this is called an atom is added
	tempAtom.Assign(element,temp.x,temp.y,temp.z);
	m_Cluster.atom.push_back(tempAtom); //push the atom on to the structure
}

/*=============================================================================
== Point_on_3D_line: takes 2 vector pairs and returns vector pair on the line
==                   defined by the 2 vector pairs nu distance along, where nu
==                   = 0 = p1 and nu = 1 = p2.
== accepts 3 points returning the third which lies between the first 2 at Nu fraction
== of the distance
=============================================================================
*/
void CGeometry::point_on_3D_line(const vec3d& b_vec, const vec3d& p2, vec3d& new_point, double nu)
{
	vec3d d_vec;

	d_vec.x = p2.x - b_vec.x;
	d_vec.y = p2.y - b_vec.y;
	d_vec.z = p2.z - b_vec.z;

	new_point.x = b_vec.x + nu * d_vec.x;
	new_point.y = b_vec.y + nu * d_vec.y;
	new_point.z = b_vec.z + nu * d_vec.z;
}
