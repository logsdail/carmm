/* Copyright Ben Curley 2006

Project:structure_generator
Filename:
Description:


*/
#include "CIco.h"
#include <iostream>
#include <cmath>
#include <cstdlib>



using namespace std;


void CDecGen::generate_vertices(int Shell_num)
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
  //the=theb;
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

void CDecGen::generate_edges( int shell_num)
{
	edge_atoms.clear();
	vec3d new_point;
	double atom_spacing;
	atom_spacing = 1.0/shell_num;

 // int index[] = {0,1,2, 0,3,4, 2,3,10 ,0,3,11, 4,5,8, 5,6,9, 6,7,10, 4,7,11};

	for(int i = 1; i < 6; i++)
	{ 
		for(int j = 1; j < shell_num; j++)
		{  
			point_on_3D_line(vertices[0], vertices[i], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
			if(i < 5)
			{
				point_on_3D_line(vertices[i], vertices[i+1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
			else
			{
				point_on_3D_line(vertices[5], vertices[1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}
	}

	for(int i = 6; i < 11; i++)
	{ 
		for(int j = 1; j < shell_num; j++)
		{  
			point_on_3D_line(vertices[11], vertices[i], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
			if(i < 10)
			{
				point_on_3D_line(vertices[i], vertices[i+1], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
			else
			{
				point_on_3D_line(vertices[10], vertices[6], new_point, atom_spacing * j);
				edge_atoms.push_back(new_point);
			}
		}
	}

	for(int i = 1; i < 6; i++)
	{
		for(int j = 1; j < shell_num; j++)
		{  
			point_on_3D_line(vertices[i], vertices[i+5], new_point, atom_spacing * j);
			edge_atoms.push_back(new_point);
		}
	}
}

void CDecGen::generate_face111( int shell_num)
{
	single_face.clear();
	edge_atoms.clear(); //this is to ensure empty vector    
	vec3d new_point;
	double atom_spacing = 1.0/shell_num; //this is more variable than in the case with edges
	int k_minus_one = shell_num - 1;


  //make 2 edges for single face
	for(int j = 1; j < shell_num; j++)
	{  
		point_on_3D_line(vertices[0], vertices[1], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

	for(int j = 1; j < shell_num; j++)
	{    
		point_on_3D_line(vertices[0], vertices[2], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

  //now build face inbetween the two edges
	for(int i = 1 ; i < shell_num-1; i++)
	{
		atom_spacing = 1.0/(i+1.0); 
		for(int j = 1; j < i+1; j++)
		{
			point_on_3D_line(edge_atoms[i], edge_atoms[i+(k_minus_one)], new_point, atom_spacing * j);
			single_face.push_back(new_point);
		}
	}

  //rotate this face 4 times in z axis
    //a 72degree rotation of these atoms in the z axis 5 times

	for(int i = 0; i < single_face.size(); i++)
	{
		face_atoms.push_back(single_face[i]);
	}
	vector<vec3d> temp, temp2; 
	temp = single_face;
	temp2 = single_face;
  //rotation matrix for each point should sort this out.
	double updateY, updateX;
	float radians = 1.25663706;
	for(int i = 1; i < 5; i++)
	{
		for(int j = 0; j < temp.size(); j++)
		{
			updateX = temp[j].x * cos(radians*i) + temp[j].y * sin(radians*i);
			updateY = temp[j].y * cos(radians*i) - temp[j].x * sin(radians*i);
			temp[j].y = updateY;
			temp[j].x = updateX;
		}

		for(int j = 0; j < temp.size(); j++)
		{
			face_atoms.push_back(temp[j]);
		}
		temp = temp2;
	}
	edge_atoms.clear();
	temp.clear();
	single_face.clear();
	atom_spacing = 1.0/shell_num; 
  //make 2 edges for single face
	for(int j = 1; j < shell_num; j++)
	{  
		point_on_3D_line(vertices[11], vertices[6], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

	for(int j = 1; j < shell_num; j++)
	{    
		point_on_3D_line(vertices[11], vertices[7], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

  //now build face inbetween the two edges
	for(int i = 1 ; i < shell_num-1; i++)
	{
		atom_spacing = 1.0/(i+1.0); 
		for(int j = 1; j < i+1; j++)
		{
			point_on_3D_line(edge_atoms[i], edge_atoms[i+(k_minus_one)], new_point, atom_spacing * j);
			single_face.push_back(new_point);
		}
	}
	temp = single_face;
	temp2 = single_face;
	for(int i = 0; i < temp.size(); i++)
	{
		face_atoms.push_back(temp[i]);
	}
	radians = 1.25663706;
	for(int i = 1; i < 5; i++)
	{
		for(int j = 0; j < temp.size(); j++)
		{
			updateX = temp[j].x * cos(radians*i) + temp[j].y * sin(radians*i);
			updateY = temp[j].y * cos(radians*i) - temp[j].x * sin(radians*i);
			temp[j].y = updateY;
			temp[j].x = updateX;
		}

		for(int j = 0; j < temp.size(); j++)
		{
			face_atoms.push_back(temp[j]);
		}
		temp = temp2;
	}
	edge_atoms.clear();
	single_face.clear();
}

void CDecGen::generate_face100( int shell_num)
{
	edge_atoms.clear();
	single_face.clear();
  //generate 1 face then rotate 4 copies is probably best approach
	double atom_spacing;
	atom_spacing = 1.0/shell_num;
	vec3d new_point;

  //this generates the edges for a single face;
	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[1], vertices[2], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[2], vertices[7], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[7], vertices[6], new_point, atom_spacing * j);
		edge_atoms.push_back(new_point);
	}

	for(int j = 1; j < shell_num; j++)
	{
		point_on_3D_line(vertices[6], vertices[1], new_point, atom_spacing * j);
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

  //single face generated
  //edge atoms are now finished so clear the edge vector 
	edge_atoms.clear();
	for(int i = 0; i < single_face.size(); i++)
	{
		face_atoms.push_back(single_face[i]);
	}

  //a 72degree rotation of these atoms in the z axis 5 times

	vector<vec3d> temp, temp2; 
	temp = single_face;

  //rotation matrix for each point should sort this out.
	double updateY, updateX;
	float radians = 1.25663706;
	for(int i = 1; i < 5; i++)
	{
		for(int j = 0; j < temp.size(); j++)
		{
			updateX = temp[j].x * cos(radians*i) + temp[j].y * sin(radians*i);
			updateY = temp[j].y * cos(radians*i) - temp[j].x * sin(radians*i);
			temp[j].y = updateY;
			temp[j].x = updateX;
		}

		for(int j = 0; j < temp.size(); j++)
		{
			face_atoms.push_back(temp[j]);
		}
		temp = single_face;
	}

}

void CDecGen::generateCoordinates(const int &shell_num, std::string element)
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
			add_coordinate(temp, element);
			generate_vertices(shell_num);
			combine_arrays(element);
			break;
		case 2:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp,element);
			generate_vertices(shell_num-1);
			combine_arrays(element);
			vertices.clear();
			generateSingleShell(shell_num,element);
			combine_arrays(element);
			break;
		case 90:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp,element);
			generate_vertices(shell_num-2);
			combine_arrays(element);
			vertices.clear();
			generateSingleShell(shell_num-1,element);
			combine_arrays(element);
			generateSingleShell(shell_num, element);
			combine_arrays(element);
			break;			
		default:
			temp.x = 0.0;
			temp.y = 0.0;
			temp.z = 0.0;
			add_coordinate(temp, element);
			generate_vertices(1);
			combine_arrays(element);
			for(int i = 2; i < shell_num+1; i++)
			{
				generateSingleShell(i,element);

			}
			break;
	}


}

void CDecGen::generateSingleShell(const int &shell_num, std::string element)
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
			generate_face111(shell_num);
			generate_edges(shell_num);
			combine_arrays(element);
	}
}
void CDecGen::combine_arrays(std::string element)
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
	vertices.clear();
	edge_atoms.clear();
	face_atoms.clear();
	single_face.clear();
};
