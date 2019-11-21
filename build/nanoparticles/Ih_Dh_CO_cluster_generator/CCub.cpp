/***************************************************************************
 *   Copyright (C) 2006 by Benjamin Curley                                 *
 *   curley@tc.bham.ac.uk						   *	
 *   rhodan@blueyonder.co.uk		                                   *
 * 									   *
 *   GeometryMaker: Creates simple 12 vectrex geometries for metallic      *
 *		    clusters						   *	                                                                    *
 ***************************************************************************/

#include "CIco.h"
#include <cstdlib>

using namespace std;

void CCubGen::generateCoordinates(const int &shell_num, std::string element)
{
	vec3d temp;
	if(shell_num < 0)
	{
		cout << "invalid Shell number must not be a negative integer, exiting" << endl;
		exit(1);
	}
	switch(shell_num)
	{
		case 0:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			break;
		case 1:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp,element);
			generate_vertices(shell_num);
			combine_arrays(element);
			break;
		case 2:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			generate_vertices(shell_num-1);
			combine_arrays(element);
			vertices.clear();
			generateSingleShell(shell_num, element);
			combine_arrays(element);
			break;
		/*		case 90: //oh what this?
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			generate_vertices(shell_num-2);
			combine_arrays(element);
			vertices.clear(element);
			generateSingleShell(shell_num-1,element);
			combine_arrays(element);
			generateSingleShell(shell_num, element);
			combine_arrays(element);
			break;	*/		
		default:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			generate_vertices(1);
			combine_arrays(element);
			for(int i = 2; i < shell_num+1; i++)
			{
				generateSingleShell(i, element);

			}
			break;
	}

}

void CCubGen::generateSingleShell(const int &shell_num, std::string element)
{
	vec3d temp;
	switch(shell_num)
	{
		case 0:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			break;
		case 1:
			generate_vertices(shell_num);
			combine_arrays(element);
			break;
		default:
			generate_vertices(shell_num);
			generate_face100(shell_num);
			edge_atoms.clear();//reset edges for 111 face creation
			generate_edges(shell_num);
			generate_face111(shell_num);
			combine_arrays(element);
	}
}

void CCubGen::generate_vertices(int shell_num)
{
	double r = shell_num * m_fRadius;
	vec3d coords;

  //Top Atom square face

	coords.x = 0.0 * r;
	coords.y = 0.5 * r;
	coords.z = 0.5 * r;
	vertices.push_back(coords);

	coords.x = 0.5 * r;
	coords.y = 0.0 * r;
	coords.z = 0.5 * r;
	vertices.push_back(coords);

	coords.x = 0.0 * r;
	coords.y =-0.5 * r;
	coords.z = 0.5 * r;
	vertices.push_back(coords);

	coords.x =-0.5 * r;
	coords.y = 0.0 * r;
	coords.z = 0.5 * r;
	vertices.push_back(coords);

  //Bottom Atom square face

	coords.x =  0.0 * r;
	coords.y =  0.5 * r;
	coords.z = -0.5 * r;
	vertices.push_back(coords);

	coords.x =  0.5 * r;
	coords.y =  0.0 * r;
	coords.z = -0.5 * r;
	vertices.push_back(coords);

	coords.x =  0.0 * r;
	coords.y = -0.5 * r;
	coords.z = -0.5 * r;
	vertices.push_back(coords);

	coords.x = -0.5 * r;
	coords.y =  0.0 * r;
	coords.z = -0.5 * r;
	vertices.push_back(coords);

  //Middle Atom square face

	coords.x = 0.5 * r;
	coords.y = 0.5 * r;
	coords.z = 0.0 * r;
	vertices.push_back(coords);

	coords.x =  0.5 * r;
	coords.y = -0.5 * r;
	coords.z =  0.0 * r;
	vertices.push_back(coords);

	coords.x = -0.5 * r;
	coords.y = -0.5 * r;
	coords.z =  0.0 * r;
	vertices.push_back(coords);

	coords.x = -0.5 * r;
	coords.y =  0.5 * r;
	coords.z =  0.0 * r;
	vertices.push_back(coords);

}

void CCubGen::generate_face100(int shell_num)
{
	double atom_spacing;
	atom_spacing = 1.0/shell_num;

	vec3d new_point;

	//this generates the edges for a single face;
	for(int i = 0; i < 3; i++ )
	{ 
		for(int j = 1; j < shell_num; j++)
		{  
			point_on_3D_line(vertices[i], vertices[i+1], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}

	for(int j = 1; j < shell_num; j++)
	{  
		point_on_3D_line(vertices[3], vertices[0], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}


	//this generates a single face
	for(int i = 0; i < shell_num-1; i++) //as we start at zero not 1 so k-2 not k-1
	{
		atom_spacing = 1.0/shell_num;
		for(int j = 1; j < shell_num; j++)
		{
			point_on_3D_line(edge_atoms[i],edge_atoms[((2*(shell_num-1))+((shell_num-2) - i))], new_point, atom_spacing * j);
			single_face.push_back(new_point);
		}

	}

	for(int i = 0; i < single_face.size(); i++)
	{
		face_atoms.push_back(single_face[i]);
	}

	//translate the face to the bottom one
	for(int i = 0; i < single_face.size(); i++)
	{
		single_face[i].z = vertices[5].z;
	}

	//now add to face_atoms vector
	for(int i = 0; i < single_face.size(); i++)
	{
		face_atoms.push_back(single_face[i]);
	}

	//a 90 degree rotation of these atoms in the Y axis face_atoms_array would lead to 2 more face complete
	//we need to make a copy

	vector<vec3d> tempY, tempX;
	tempY = face_atoms;
	tempX = face_atoms;

	//rotation matrix for each point should sort this out.

	//rotation in the x plane for tempX
	double updateY, updateZ;
	for(int i = 0; i < tempX.size(); i++)
	{
		updateY = tempX[i].y * cos(1.57079633) + tempX[i].z * sin(1.57079633); 
		updateZ = tempX[i].z * cos(1.57079633) - tempX[i].y * sin(1.57079633);
		tempX[i].y = updateY;
		tempX[i].z = updateZ;

	}

	//rotation in y plane for tempY 1.5709663 = 90 degrees in radians
	double updateX;

	for(int i = 0; i < tempY.size(); i++)
	{
		updateX = tempY[i].x * cos(1.57079633) - tempY[i].z * sin(1.57079633); 
		updateZ = tempY[i].x * sin(1.57079633) + tempY[i].z * cos(1.57079633);
		tempY[i].x = updateX;
		tempY[i].z = updateZ;
	}

	for(int i = 0; i < tempX.size(); i++)
	{
		face_atoms.push_back(tempX[i]);
	}

	for(int i = 0; i < tempY.size(); i++)
	{
		face_atoms.push_back(tempY[i]);
	}
}


void CCubGen::generate_face111(int shell_num)
{

	vec3d new_point;
	double atom_spacing; //this is more variable than in the case with edges
	int k_minus_one = shell_num - 1;

  //8 of the 8 111 faces
	int step = 3 * (shell_num-1); 
	for(int k = 0; k < 8; k++)
	{ 
		for(int i = 3 ; i <  step; i+= 3)
		{
			atom_spacing = 1.0/((i/3.0)+1.0); //some what over complicated but it brings the 3 back into integers and adds the approiatte amount for correct divison.
			int number_required = (i/3)+1;
			for(int j = 1; j < number_required; j++)
			{
				point_on_3D_line(edge_atoms[i + ( step * k)], edge_atoms[i+1 +(step * k)], new_point, atom_spacing * j);
				face_atoms.push_back(new_point);
			}
		}
	}

}

void CCubGen::generate_edges(int shell_num)
{
	vec3d new_point;
	double atom_spacing;
	atom_spacing = 1.0/shell_num;

	int index[] = {0,1,8, 1,2,9, 2,3,10 ,0,3,11, 4,5,8, 5,6,9, 6,7,10, 4,7,11};

	for(int i = 0; i < 22; i += 3)
	{ 
		for(int j = 1; j < shell_num; j++)
		{  
			point_on_3D_line(vertices[index[i]], vertices[index[i+1]], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);

			point_on_3D_line(vertices[index[i]], vertices[index[i+2]], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);

			point_on_3D_line(vertices[index[i+1]], vertices[index[i+2]], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}
}

void CCubGen::combine_arrays(std::string element)
{
	vec3d temp;
	for(int i = 0; i < vertices.size(); i++)
	{
		temp = vertices[i];
		add_coordinate(temp, element);
	}
	for(int i = 0; i < edge_atoms.size(); i++)
	{
		temp = edge_atoms[i];
		add_coordinate(temp, element);
	}
	for(int i = 0; i < face_atoms.size(); i++)
	{
		temp = face_atoms[i];
		add_coordinate(temp, element);
	}

	vertices.clear();
	edge_atoms.clear();
	face_atoms.clear();
	single_face.clear();

}
