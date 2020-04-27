/***************************************************************************
 *   Copyright (C) 2005 by Ben Curley                                      *
 *   curley@tc.bham.ac.uk                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#pragma once;
#include <sstream>
#include <iostream>
#include <fstream> //stl file handling
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

struct atom
{
	string atom_type;
	double x;
	double y;
	double z;
	int coordination;
	int identity;
};

struct plane
{
	atom *patom_1;
	atom *patom_2;
	atom *patom_3;
	bool surface;
};

struct atomData
{
	std::string atomType;
	double x;
	double y;
	double z;
	void Assign(std::string name, double X, double Y, double Z)
	{
		atomType = name;
		x = X;
		y = Y;
		z = Z;
	}
};

struct clusterData
{
	int size; //Number of atoms
	std::string comment; //commentline
	double potentialEnergy; //potential energy of cluster
	std::vector<atomData> atom;

	clusterData() : size(0), comment(""), potentialEnergy(0)
	{
	};	
	friend std::ostream& operator<<(std::ostream& os, const clusterData& c);
};
