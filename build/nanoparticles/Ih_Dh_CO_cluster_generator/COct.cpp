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

void COctGen::generateCoordinates(const int &shell_num, std::string element)
{
	if(shell_num < 0)
	{
		cout << "invalid Shell number must not be a negative integer, exiting" << endl;
		exit(1);
	}
	switch(shell_num)
	{
		case 0:
			generateOrigin(shell_num,element);
			break;
		case 1:
			generateOrigin(shell_num,element);
			generate_vertices(shell_num);
			combine_arrays(shell_num,element);
			break;
		case 2:
			generateOrigin(shell_num,element);
			generate_vertices(shell_num-1);
			combine_arrays(shell_num-1,element);
			vertices.clear();
			generateSingleShell(shell_num, element);
			combine_arrays(shell_num,element);
			break;
		default:
			generateOrigin(shell_num,element);
			generate_vertices(1);
			combine_arrays(shell_num,element);
			for(int i = 2; i < shell_num+1; i++)
			{
				generateSingleShell(i, element);

			}
			break;
	}

}

void COctGen::generateSingleShell(const int &shell_num, std::string element)
{
	switch(shell_num)
	{
		case 0:
			generateOrigin(shell_num,element);
			break;
		case 1:
			generate_vertices(shell_num);
			combine_arrays(shell_num,element);
			break;
		default:
			generate_vertices(shell_num);
			generate_edges(shell_num);
			generate_face100(shell_num);
			//generate_edges(shell_num);
			combine_arrays(shell_num,element);
	}
}

void COctGen::generateOrigin(int shell_num, std::string element)
{
	vec3d temp;
	temp.x = 0.0;
        temp.y = 0.0;
        temp.z = 0.0;
        add_coordinate(temp, element);

	if (shell_num != m_iShells)
        {
        	temp.z = m_fRadius * m_iShells;
                add_coordinate(temp, element);
        }
}

void COctGen::generate_vertices(int shell_num)
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
}

void COctGen::generate_face100(int shell_num)
{
	double atom_spacing;
	atom_spacing = 1.0/shell_num;

	vec3d new_point;

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
}

void COctGen::generate_edges(int shell_num)
{
	vec3d new_point;
        double atom_spacing;
        atom_spacing = 1.0/shell_num;

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
}

void COctGen::combine_arrays(int shell_num, std::string element)
{
	vec3d temp;

	for(int i = 0; i < vertices.size(); i++)
	{
		temp = vertices[i];
		add_coordinate(temp, element);

		if (shell_num != m_iShells)
		{
			temp.z = 0.5 * m_fRadius * ((2 * m_iShells) - shell_num);
			add_coordinate(temp, element);
		}
        } 

	for(int i = 0; i < edge_atoms.size(); i++)
	{
		temp = edge_atoms[i];
		add_coordinate(temp, element);

                if (shell_num != m_iShells)
                {
                        temp.z = 0.5 * m_fRadius * ((2 * m_iShells) - shell_num);
                        add_coordinate(temp, element);
                }
	}
	for(int i = 0; i < face_atoms.size(); i++)
	{
		temp = face_atoms[i];
		add_coordinate(temp, element);

                if (shell_num != m_iShells)
                {
                        temp.z = 0.5 * m_fRadius * ((2 * m_iShells) - shell_num);
                        add_coordinate(temp, element);
                }
	}

	vertices.clear();
	edge_atoms.clear();
	face_atoms.clear();
	single_face.clear();
}
