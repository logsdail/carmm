/***************************************************************************
 *   Copyright (C) 2006 by Benjamin Curley                                 *
 *   curley@tc.bham.ac.uk						  						   *	
 *   rhodan@blueyonder.co.uk		                                   	   *
 * 									                                       *
 *   GeometryMaker: Creates simple 12 vectrex geometries for metallic      *
 *		    clusters						     	                       *
 ***************************************************************************/

 #include "CIco.h"
 #include <iostream>
 #include <cmath>



using namespace std;

CIcoGen::~CIcoGen()
{

}

CIcoGen::CIcoGen()
{
}

void CIcoGen::generateCoordinates(const int &shell_num, std::string element)
{
	vec3d point;
	if( shell_num == 0 )
	{
		point.x = 0.0;
		point.y = 0.0;
		point.z = 0.0;
		add_coordinate(point, element);
		clear_arrays();
	}
	if( shell_num == 1 )
	{
		point.x = 0.0;
		point.y = 0.0;
		point.z = 0.0;
		add_coordinate(point,element);
		clear_arrays();
		generate_vertices(shell_num);
		combine_arrays(element);
		clear_arrays();
	}
	if( shell_num == 2 )
	{
		point.x = 0.0;
		point.y = 0.0;
		point.z = 0.0;
		add_coordinate(point,element);
		clear_arrays();
		generate_vertices(1);
		combine_arrays(element);
		clear_arrays();
		generate_vertices(shell_num);
		generate_edges(shell_num);
		combine_arrays(element);
		clear_arrays();
	}
	if( shell_num > 2 )
	{
		point.x = 0.0;
		point.y = 0.0;
		point.z = 0.0;
		add_coordinate(point,element);
		clear_arrays();
		generate_vertices(1);
		combine_arrays(element);
		clear_arrays();
		generate_vertices(2);
		generate_edges(2);
		combine_arrays(element);
		clear_arrays();
		for(int i = 3; i < shell_num+1; i++)
		{
			generate_vertices(i);
			generate_edges(i);
			generate_faces(i);
			combine_arrays(element);
			clear_arrays();
		}
	}

}


void CIcoGen::generateSingleShell(const int &shell_num, std::string element)
{
	vec3d point;
	if( shell_num == 0 )
	{
		point.x = 0.0;
		point.y = 0.0;
		point.z = 0.0;
		add_coordinate(point,element);
		clear_arrays();
	}
	if( shell_num == 1 )
	{
		generate_vertices(shell_num);
		combine_arrays(element);
		clear_arrays();
	}
	if( shell_num == 2 )
	{
		generate_vertices(shell_num);
		generate_edges(shell_num);
		combine_arrays(element);
		clear_arrays();
	}
	if( shell_num > 2 )
	{
		generate_vertices(shell_num);
		generate_edges(shell_num);
		generate_faces(shell_num);
		combine_arrays(element);
		clear_arrays();
	}

}


void CIcoGen::combine_arrays(std::string element)
{
	vec3d temp;
	for(int i = 0; i < vertices.size(); i++)
	{
		temp = vertices[i];
		add_coordinate(temp,element);
	}
	for(int i = 0; i < edge_atoms.size(); i++)
	{
		temp = edge_atoms[i];
		add_coordinate(temp,element);
	}
	for(int i = 0; i < face_atoms.size(); i++)
	{
		temp = face_atoms[i];
		add_coordinate(temp,element);
	}

}

/*
========================================================================
== generate_ico_vertices: places the 12 points of and Icosahedron within
==                        a sphere of radius r.
========================================================================
*/
void CIcoGen::generate_vertices(int Shell_num)
{
	double r = Shell_num * m_fRadius;
	vec3d coords;
	double phiaa = 26.5605;
	double phia = PI*phiaa/180.0;
	double theb = PI*36.0/180.0;
	double the72 = PI*72.0/180.0;

  //Top Atom

	coords.x = 0.0;
	coords.y = 0.0;
	coords.z = r;
	vertices.push_back(coords);

	double the = 0.0;
	for(int i=1; i < 6; i++)
	{
		coords.x = r * cos(the) * cos(phia);
		coords.y = r * sin(the) * cos(phia);
		coords.z = r * sin(phia);
		the += the72;
		vertices.push_back(coords);
	}
	the=theb;
	for(int i=6; i < 11; i++)
	{
		coords.x = r * cos(the) * cos(-phia);
		coords.y = r * sin(the) * cos(-phia);
		coords.z = r * sin(-phia);
		the += the72;
		vertices.push_back(coords);
	}

  //Bottom Atom
	coords.x = 0.0;
	coords.y = 0.0;
	coords.z = -r;
	vertices.push_back(coords);
	if( vertices.size() != 12)
	{
		cout << "Shit" << endl;
	}
}

void CIcoGen::generate_faces(int shell_num)
{
	if( shell_num < 3 )
	{
		cout << "The cluster shell is to small for face atoms to exist" << endl;
		return;
	}

	vec3d new_point;
	double atom_spacing; //this is more variable than in the case with edges
	int k_minus_one = shell_num - 1;
  // 1 < 3
	for(int j = 1; j < k_minus_one; j++) //this moves down the edge atoms shell_num -1 = number of edge atoms on the edge, the first atom on an edge is to close
                                        //to its neighbouring atom for the space so should start at 1?
	{
		atom_spacing = 1.0/(j+1);
		for(int i = 0; i < 5; i++) //this loop moves around the 5 faces "well 4 atm"
		{                           //1+(0*3)/1+(1*3)/1+(2*3)    //1+(1*3)/1+(2*3)
    //cout << j+(i*(shell_num-1)) << endl;
			for(int k = 1; k < j+1; k++)
			{
				if( i != 4 )
				{
					point_on_3D_line(edge_atoms[j + i * k_minus_one], edge_atoms[j+ (i+1) * k_minus_one], new_point, atom_spacing * k);
					face_atoms.push_back(new_point);
				}
				else
				{
					point_on_3D_line(edge_atoms[j + i * k_minus_one], edge_atoms[j + 0], new_point, atom_spacing * k);
					face_atoms.push_back(new_point);
				}
			}
		}
	}

  //Repeat for the lower Half of the Icosahredron but starting at the right place (5*k-1 should be the right place)
	for(int j = 1; j < k_minus_one; j++) //this moves down the edge atoms shell_num -1 = number of edge atoms on the edge, the first atom on an  edge is to close
                                        //to its neighbouring atom for the space so should start at 1?
	{
		atom_spacing = 1.0/(j+1);
		for(int i = 0; i < 5; i++) //this loop moves around the 5 faces "well 4 atm"
		{                           //1+(0*3)/1+(1*3)/1+(2*3)    //1+(1*3)/1+(2*3)
    //cout << j+(i*(shell_num-1)) << endl;
			for(int k = 1; k < j+1; k++)
			{
				if( i != 4 )
				{
					point_on_3D_line(edge_atoms[j + i * k_minus_one +  5 * k_minus_one], 
						edge_atoms[j + (i+1) * k_minus_one + 5 * k_minus_one], 
						new_point, atom_spacing * k
					);
					face_atoms.push_back(new_point);
				}
				else
				{
					point_on_3D_line(edge_atoms[j + i * k_minus_one + 5 * k_minus_one], 
						edge_atoms[j + 0 + 5 * k_minus_one], 
						new_point, atom_spacing * k
					);
					face_atoms.push_back(new_point);
				}
			}
		}
	}
  //we start on (k-1)*10
  //

	for(int h = 0; h < k_minus_one; h++)
	{
		atom_spacing = 1.0/(k_minus_one-h);
		for(int i = 0; i < 4; i++) //for 4 faces the 5th face special case
		{
			for(int j = 1; j < (k_minus_one-h)/*thisneeds generalising*/; j++)
			{
				point_on_3D_line(edge_atoms[(k_minus_one * 10 + h) + (k_minus_one * i)], edge_atoms[(k_minus_one * 15 + k_minus_one + h) 
					+ (k_minus_one * i)], new_point, atom_spacing * j);
				face_atoms.push_back(new_point);
			}
		}
	}
    //the final face for \/ faces
	for(int h = 0; h < k_minus_one; h++)
	{
		atom_spacing = 1.0/(k_minus_one-h);
		for(int i = 4; i < 5; i++) //for 4 faces the 5th face special case
		{
			for(int j = 1; j < (k_minus_one-h)/*thisneeds generalising*/; j++)
			{//this should = 3*10+0+3*4 = 42 which it is and wants to get to 45 for k=4 case 3*15+3+0 = 48
				point_on_3D_line(edge_atoms[(k_minus_one * 10 + h) + (k_minus_one * i)], 
					edge_atoms[(k_minus_one * 15 + h)], 
					new_point, atom_spacing * j);
				face_atoms.push_back(new_point);
			}
		}
	}
// Now for the faces /\ in the icosahedron
	for(int h = 0; h < k_minus_one; h++)
	{
		atom_spacing = 1.0/(k_minus_one-h);
		for(int i = 0; i < 4; i++) //for 4 faces the 5th face special case
		{
			for(int j = 1; j < (k_minus_one-h)/*thisneeds generalising*/; j++)
			{// 45 + 2 + 0 = 47 45 + 2 + 3 = 50 45 + 2 + 6 = 53 //45+2+3-1 = 49
				point_on_3D_line(edge_atoms[(k_minus_one * 15 + (k_minus_one-1)) + (k_minus_one * i) - h],
					edge_atoms[(k_minus_one * 10 + (k_minus_one-1) - h) + (k_minus_one * i)], 
					new_point, atom_spacing * j
				);
				face_atoms.push_back(new_point);
			}
		}
	}

	for(int h = 0; h < k_minus_one; h++)
	{
		atom_spacing = 1.0/(k_minus_one-h);
		for(int i = 4; i < 5; i++) //for 4 faces the 5th face special case
		{
			for(int j = 1; j < (k_minus_one-h)/*thisneeds generalising*/; j++)
			{// 45+2+12-0=59
				point_on_3D_line(edge_atoms[(k_minus_one * 15 + (k_minus_one-1)) + (k_minus_one * i) - h], edge_atoms[((k_minus_one) * 15 - 1
						- h)], new_point, atom_spacing * j);
				face_atoms.push_back(new_point);
			}
		}
	}

}

void CIcoGen::generate_edges(int shell_num)
{
  //atoms on edge are related by shellsize k-1 = num atoms per edges
  //k=3 then atoms = 2 vert-----at-------at------vert so the edge length should
  //be divded by shell_size 1.0/k * j where j >= 1 <= k-1
	vec3d new_point;
	double atom_spacing;
	atom_spacing = 1.0/shell_num;

  //Place atoms along edges that connect top of the ico to the 1st pentagon
	for(int i = 1; i < 6; i++)
	{
		for(int j = 1; j < shell_num; j++)
		{
			point_on_3D_line(vertices[0], vertices[i], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}
  //Number of atoms currently occupying the edge_atoms vector = (k-1) *5 "thats for 5 edges"
  //cout << edge_atoms.size() << " = " << (shell_num-1)*5 << endl;
  //Place atoms along edges that connect bottom of the ico to the 2nd pentagon
	for(int i = 6; i < 11; i++)
	{
		for(int j = 1; j < shell_num; j++)
		{
			point_on_3D_line(vertices[11], vertices[i], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}

  // The array should now be of size (k-1) * 10
  //cout << edge_atoms.size() << " = " << (shell_num-1)*10 << endl;
  // put atoms along the lines connecting each pentagon, misses the final one

	for(int i = 2; i < 6; i++)
	{
		for(int j = 1; j < shell_num; j++)
		{
			point_on_3D_line(vertices[i], vertices[i+5], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}

	}
	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[1], vertices[6], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

    //point_on_3D_line(vertices[i], vertices[i+4], new_point, atom_spacing * j);
    //edge_atoms.push_back(new_point);
    // }
    //}

    //the array should now be the size of k-1 * 15

	for(int i = 2; i < 6; i++)
	{
		for(int j = 1; j < shell_num; j++)
		{
			point_on_3D_line(vertices[i], vertices[i+4], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}
   //this is the final one
	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[1], vertices[10], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

  //array should be k-1 * 20 in size now only 10 edges remaining
  //cout << edge_atoms.size() << " this is after we done all interconnecting triangles" << endl;

  //put atoms on edges around each pentagon
	for(int i = 1; i < 6; i++)
	{
		if( i != 5 )
		{
			for(int j = 1; j < shell_num; j++)
			{
				point_on_3D_line(vertices[i], vertices[i+1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}
		else
		{
			for(int j = 1; j < shell_num; j++)
			{
				point_on_3D_line(vertices[5], vertices[1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}
	}

  //same as above but for the lower pentagon
	for(int i = 6; i < 11; i++)
	{
		if( i != 10 )
		{
			for(int j = 1; j < shell_num; j++)
			{
				point_on_3D_line(vertices[i], vertices[i+1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}
		else
		{
			for(int j = 1; j < shell_num; j++)
			{
				point_on_3D_line(vertices[i], vertices[6], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}

	}
}
