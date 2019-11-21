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
#include "cmorse_minimiser.h"

CMorse_minimiser::CMorse_minimiser(string potential_filename)
{
	char buffer[256];
	fstream potential_file;

	potential_file.open(potential_filename.c_str(), ios::in);
	if(!potential_file.is_open())
	{
		cout << "Potential Doesnt exist" << endl;
		exit(1);
	}
	else
	{
		potential_file >> alpha;
		potential_file.getline(buffer, 256);
		potential_file >> element;
		potential_file >> bond_scale_factor;
		potential_file.getline(buffer, 256);
	}
}

CMorse_minimiser::~CMorse_minimiser()
{
}

clusterData CMorse_minimiser::GetMinimisedData()
{
	//must prepare for minimisation
	//first convert clusterData to a vector<double>
	m_inumVariables = m_Cluster.size * 3;
	m_vVariables.resize(m_inumVariables); //resize vector
	for(size_t i = 0; i < m_Cluster.size; i++)
	{
		//initalise m_vVariables
		m_vVariables[i] = m_Cluster.atom[i].x;
		m_vVariables[i + m_Cluster.size] = m_Cluster.atom[i].y;
		m_vVariables[i + m_Cluster.size * 2] = m_Cluster.atom[i].z;
	}
	//run the minimiser shall make this protected so it can't be called directly
	this->minimise();
	for(size_t i = 0; i < m_Cluster.size; i++)
	{
		m_Cluster.atom[i].x = m_vVariables[i];
		m_Cluster.atom[i].y = m_vVariables[i + m_Cluster.size];
		m_Cluster.atom[i].z = m_vVariables[i + m_Cluster.size * 2];
	}
	m_Cluster.potentialEnergy = m_dEnergy;
	return m_Cluster;
}

double CMorse_minimiser::fx_dx(const double* x, double* g)
{
	int n_atoms = m_Cluster.size;
	double dxij = 0.0;
	double dyij = 0.0;
	double dzij = 0.0;
	double rij = 0.0;

	double afac = 0.0;
	double function_eval = 0.0;
	int counter = 0;

	double dv2dr= 0.0;
	double dvdx = 0.0;
	double dvdy = 0.0;
	double dvdz = 0.0;

	int nmax = n_atoms /** 3*/ + 1; //Got to do something about this nmax

	double dv2dx[nmax];
	double dv2dy[nmax];
	double dv2dz[nmax];
	double aa = alpha*2.0;

	for(int i = 0; i < n_atoms; i++)
	{
		dv2dx[i] = 0.0;
		dv2dy[i] = 0.0;
		dv2dz[i] = 0.0;
	}



	for(int i = 0; i < n_atoms-1; i++)
	{
		for(int j = i+1; j < n_atoms; j++)
		{
			dxij = x[i] - x[j];
			dyij = x[n_atoms+i] - x[n_atoms+j];
			dzij = x[n_atoms*2+i] - x[n_atoms*2+j];
			rij = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
			afac=exp(alpha*(1.0-rij));
			function_eval += afac*(afac-2.0);

			dv2dr = aa*afac*(1.0-afac);

			double one_over_rij = 1.0 / rij;

			dvdx = dv2dr*dxij*one_over_rij;
			dvdy = dv2dr*dyij*one_over_rij;
			dvdz = dv2dr*dzij*one_over_rij;

			dv2dx[i] += dvdx;

			dv2dy[i] += dvdy;
			dv2dz[i] += dvdz;

			dv2dx[j] -= dvdx;
			dv2dy[j] -= dvdy;
			dv2dz[j] -= dvdz;
		}
	}

	for(int i = 0; i < n_atoms; i++)
	{
		g[i] = dv2dx[i];
		g[n_atoms+i] = dv2dy[i];
		g[2*n_atoms+i] = dv2dz[i];
	}

	return function_eval;

}

