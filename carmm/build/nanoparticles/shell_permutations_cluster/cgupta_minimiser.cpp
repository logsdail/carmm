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
#include "cgupta_minimiser.h"
#include "Struct.h"
#include "cError.h"

CGupta_minimiser::CGupta_minimiser(string potential_filename, bool bCutoff)
{
	m_bCutoff = bCutoff; // Defines if the Cutoff is implemented

	m_potentialFilename = potential_filename;
	this->ReadPotential();

	if (m_bCutoff){
        	this->ReadCutoff();
        }

}

CGupta_minimiser::~CGupta_minimiser()
{

}

clusterData CGupta_minimiser::GetMinimisedData()
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

void CGupta_minimiser::SetClusterData(clusterData cluster)
{
	m_Cluster = cluster;
}

//deprecated must make clear that the returned structure is the minimised one
clusterData CGupta_minimiser::GetClusterData()
{
	return m_Cluster;
}

double CGupta_minimiser::GetEnergy()
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

	this->CalculateEnergy();

	return m_dEnergy;
}

//This is a little tricker as I need to know atom types for bimetallics computations
//Implementation of potential derivatives 
double CGupta_minimiser::fx_dx(const double* x, double* g)
{
	//way to much data placed on the stack for large system so we have stack overflow
	//need to use the heap
	int n_atoms = m_Cluster.size;
	int nmax = m_inumVariables + 1; 
	double function_eval = 0.0;

	double dxij = 0.0;
	double dyij = 0.0;
	double dzij = 0.0;
	double rij = 0.0;

	double v1i = 0.0;
	double v2i = 0.0;
	double v1ij = 0.0;
	double v2ij = 0.0;
	double v2i1 = 0.0;
	double dv1dxi = 0.0;
	double dv1dyi = 0.0;
	double dv1dzi = 0.0;

	double dv2dxi = 0.0;
	double dv2dyi = 0.0;
	double dv2dzi = 0.0;

	double dv1dxij = 0.0;
	double dv1dyij = 0.0;
	double dv1dzij = 0.0;


	//shall place all of these on to the heap
	//so pointers away

	double *v1 = new double[nmax];
	double *v2 = new double[nmax];
	double *v2n = new double[nmax];
	double *dv1dx = new double[nmax];  
	double *dv1dy = new double[nmax];  
	double *dv1dz = new double[nmax];  
	double *dv2dx = new double[nmax];  
	double *dv2dy = new double[nmax];  
	double *dv2dz = new double[nmax];

	int nmax_sqr_div = int((nmax*nmax)*0.5);

	double *dv2dxij = new double[nmax_sqr_div];
	double *dv2dyij = new double[nmax_sqr_div];
	double *dv2dzij = new double[nmax_sqr_div];

	double r0ij = 0.0;
	double rhoij= 0.0;
	double aaij = 0.0;
	double z2ij = 0.0;
	double qqij = 0.0;
	double cc1ij = 0.0;
	double cc2ij = 0.0;

        double cutoff_start = 0.0;
	double cutoff_end = 0.0;
	double a3 = 0.0;
	double a4 = 0.0;
	double a5 = 0.0;
	double x3 = 0.0;
	double x4 = 0.0;
	double x5 = 0.0; 

	double rdiv = 0.0;
	double one_over_rij = 0.0;

	double AA_r0ij = 1.0/AA_r0;
	double BB_r0ij = 1.0/BB_r0;
	double AB_r0ij = 1.0/AB_r0;
	double AA_z2ij = AA_zeta * AA_zeta;
	double BB_z2ij = BB_zeta * BB_zeta;
	double AB_z2ij = AB_zeta * AB_zeta;

 // partial_energy = vector<double>(n_atoms);

	for(int i = 0; i < n_atoms; i++)
	{
		v1[i] = 0.0;
		v2[i] = 0.0;
		dv1dx[i] = 0.0;
		dv1dy[i] = 0.0;
		dv1dz[i] = 0.0;
		dv2dx[i] = 0.0;
		dv2dy[i] = 0.0;
		dv2dz[i] = 0.0;
	}

	for(int i = 0; i < nmax_sqr_div; i++)
	{
		dv2dxij[i] = 0.0;
		dv2dyij[i] = 0.0;
		dv2dzij[i] = 0.0;
	}

	int nn = 0;

	for(int i = 0; i < n_atoms; i++) 
	{
		v1i = v1[i];
		v2i = v2[i];
		dv1dxi = dv1dx[i];
		dv1dyi = dv1dy[i];
		dv1dzi = dv1dz[i];
		dv2dxi = dv2dx[i];
		dv2dyi = dv2dy[i];
		dv2dzi = dv2dz[i];

		for(int j = i+1; j < n_atoms; j++)
		{
			nn++;

			if(m_Cluster.atom[i].atomType == element_a && m_Cluster.atom[j].atomType == element_a)
			{
				r0ij = AA_r0ij;
				rhoij = AA_p;
				aaij = AA_A;
				z2ij = AA_z2ij;
				qqij = AA_q;
				cc1ij = -AA_p * r0ij;
				cc2ij = AA_q * r0ij;

        			if (m_bCutoff){
					cutoff_start = AA_cutoff_start;
					cutoff_end = AA_cutoff_end;
					a3 = AA_a3;
					a4 = AA_a4;
					a5 = AA_a5;
					x3 = AA_x3;
					x4 = AA_x4;
					x5 = AA_x5;
				}
			} 
			else if(m_Cluster.atom[i].atomType == element_b && m_Cluster.atom[j].atomType == element_b)
			{
				r0ij = BB_r0ij;
				rhoij = BB_p;
				aaij = BB_A;
				z2ij = BB_z2ij;
				qqij = BB_q;
				cc1ij = -BB_p * r0ij;
				cc2ij = BB_q * r0ij;

                              if (m_bCutoff){
                                        cutoff_start = AB_cutoff_start;
                                        cutoff_end = AB_cutoff_end;
                                        a3 = AB_a3;
                                        a4 = AB_a4;
                                        a5 = AB_a5;
                                        x3 = AB_x3;
                                        x4 = AB_x4;
                                        x5 = AB_x5;
                                }

			}
			else
			{
				r0ij = AB_r0ij;
				rhoij = AB_p;
				aaij = AB_A;
				z2ij = AB_z2ij;
				qqij = AB_q;
				cc1ij = -AB_p *  r0ij;
				cc2ij = AB_q * r0ij;

                                if (m_bCutoff){
                                        cutoff_start = BB_cutoff_start;
                                        cutoff_end = BB_cutoff_end;
                                        a3 = BB_a3;
                                        a4 = BB_a4;
                                        a5 = BB_a5;
                                        x3 = BB_x3;
                                        x4 = BB_x4;
                                        x5 = BB_x5;
                                }

			}

			dxij = x[i] - x[j];
			dyij = x[i+n_atoms] - x[j+n_atoms];
			dzij = x[i+n_atoms*2] - x[j+n_atoms*2];
			rij = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
                        
			rdiv = rij*r0ij - 1.0;
                        one_over_rij = 1.0/rij;

			double attractive = 0;
			double repulsive = 0;

			if (!m_bCutoff || rij < cutoff_start) 
			{

				v1ij = aaij * exp(-rhoij * rdiv);
				v2ij = z2ij * exp(-2.0 * qqij * rdiv);

				repulsive = v1ij * cc1ij;
				attractive = v2ij * cc2ij;

			}
			else if (rij < cutoff_end) 
			{

	                        double factor1, factor2, factor3, factor4, factor5;

                                factor1 = rij-cutoff_end;
                                factor2 = factor1*factor1;
                                factor3 = factor2*factor1;
                                factor4 = factor3*factor1;
                                factor5 = factor4*factor1;

                                v1ij = a5*factor5 + a4*factor4 + a3*factor3;
                                v2ij = x5*factor5 + x4*factor4 + x3*factor3;

				attractive = - ((5*x5*factor4 + 4*x4*factor3 + 3*x3*factor2)*v2ij);
				repulsive = (5*a5*factor4 + 4*a4*factor3 + 3*a3*factor2);

				v2ij = v2ij*v2ij;

			} 
			else
			{
				
                                v1ij = 0;
                                v2ij = 0;

				repulsive = 0;
				attractive = 0;

			}

	// End of current Andrew Logsdail additions....

			v1i += v1ij;
			v2i += v2ij;

			v1[j] += v1ij;
			v2[j] += v2ij;

	  //Deravitives of V1

                        dv1dxij = repulsive * dxij * one_over_rij;
                        dv1dyij = repulsive * dyij * one_over_rij;
                        dv1dzij = repulsive * dzij * one_over_rij;

			dv1dxi += dv1dxij;
			dv1dyi += dv1dyij;
			dv1dzi += dv1dzij;

			dv1dx[j] -= dv1dxij;
			dv1dy[j] -= dv1dyij;
			dv1dz[j] -= dv1dzij;

	  //derivatives of V2

                        dv2dxij[nn] = attractive * dxij * one_over_rij;
                        dv2dyij[nn] = attractive * dyij * one_over_rij;
                        dv2dzij[nn] = attractive * dzij * one_over_rij;

			dv2dxi += dv2dxij[nn];
			dv2dyi += dv2dyij[nn];
			dv2dzi += dv2dzij[nn];

			dv2dx[j] -= dv2dxij[nn];
			dv2dy[j] -= dv2dyij[nn];
			dv2dz[j] -= dv2dzij[nn];

		}

		v2i1 = 1.0/sqrt(v2i);


		function_eval += (v1i - sqrt(v2i));
      //partial_energy[i] = (v1i -sqrt(v2i));

		v2n[i] = v2i1; //was v2[i] but this was becuase its not used outside j loop but I dont lie that

		dv1dx[i] = dv1dxi;
		dv1dy[i] = dv1dyi;
		dv1dz[i] = dv1dzi;

		dv2dx[i] = v2i1 * dv2dxi;
		dv2dy[i] = v2i1 * dv2dyi;
		dv2dz[i] = v2i1 * dv2dzi;
	}

	nn = 0;

	for(int i = 0; i < n_atoms; i++)
	{
		for(int j = i+1; j < n_atoms; j++)
		{
			nn++;
			dv2dx[j] -= v2n[i] * dv2dxij[nn];
			dv2dy[j] -= v2n[i] * dv2dyij[nn];
			dv2dz[j] -= v2n[i] * dv2dzij[nn];

			dv2dx[i] += v2n[j] * dv2dxij[nn];
			dv2dy[i] += v2n[j] * dv2dyij[nn];
			dv2dz[i] += v2n[j] * dv2dzij[nn];
		}
	}

	for(int i = 0; i < n_atoms; i++)
	{
		g[i] = 2.0 * dv1dx[i] + dv2dx[i];
		g[i+n_atoms] = 2.0 * dv1dy[i] + dv2dy[i];
		g[i+2*n_atoms] = 2.0 * dv1dz[i] + dv2dz[i];
	} 

	//ok we dont everything best delete the variables on the stack

	delete [] v1; 
	delete [] v2;
	delete [] v2n;
	delete [] dv1dx;  
	delete [] dv1dy;  
	delete [] dv1dz;  
	delete [] dv2dx;  
	delete [] dv2dy;  
	delete [] dv2dz;

	delete [] dv2dxij;
	delete [] dv2dyij;
	delete [] dv2dzij;

//	This was used for testing when inserting cutoff function
//	std::cout << function_eval << std::endl;

	return function_eval;

}

void CGupta_minimiser::ReadCutoff()
{
        std::fstream cutoff;
        std::string inputname = "cutoff_parameters";
        char buffer[256];

	cutoff.open(inputname.c_str(), std::ios::in);
        if(!cutoff.is_open())
        {
                throw cError ("Unable to open cutoff_parameters!");
        }
	else
	{
                
                cutoff >> num_els_cutoff;

		if ( num_els_cutoff != num_els )
		{
			throw cError ("Number of elements in potential and cutoff do not match");
		}

		// get cutoff parameters
                cutoff >> AA_cutoff_start;
                cutoff.getline(buffer, 256);
                cutoff >> AA_cutoff_end;
                cutoff.getline(buffer, 256);
                cutoff >> AA_a3;
                cutoff.getline(buffer, 256);
                cutoff >> AA_a4;
                cutoff.getline(buffer, 256);
                cutoff >> AA_a5;
                cutoff.getline(buffer, 256);
                cutoff >> AA_x3;
                cutoff.getline(buffer, 256);
                cutoff >> AA_x4;
                cutoff.getline(buffer, 256);
                cutoff >> AA_x5;
                cutoff.getline(buffer, 256);

		if ( num_els_cutoff == 2 )
		{
                        cutoff >> AB_cutoff_start;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_cutoff_end;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_a3;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_a4;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_a5;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_x3;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_x4;
                        cutoff.getline(buffer, 256);
                        cutoff >> AB_x5;
                        cutoff.getline(buffer, 256);

			cutoff >> BB_cutoff_start;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_cutoff_end;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_a3;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_a4;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_a5;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_x3;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_x4;
                	cutoff.getline(buffer, 256);
                	cutoff >> BB_x5;
                	cutoff.getline(buffer, 256);
		}

	}

	cutoff.close();

}


void CGupta_minimiser::ReadPotential()
{
	fstream potential_file;
	char buffer[256];

	potential_file.open(m_potentialFilename.c_str(), ios::in);
	if(!potential_file.is_open())
	{
		throw cError(" Unable to open potential file!");
	}
	else
	{
		potential_file >> num_els;
		potential_file >> AA_A;
		potential_file.getline(buffer, 256);
		potential_file >> AA_zeta;
		potential_file.getline(buffer, 256);
		potential_file >> AA_p;
		potential_file.getline(buffer, 256);
		potential_file >> AA_q;
		potential_file.getline(buffer, 256);
		potential_file >> AA_r0;
		potential_file.getline(buffer, 256);

		if ( num_els == 2 )
		{
			potential_file >> AB_A;
			potential_file.getline(buffer, 256);
			potential_file >> AB_zeta;
			potential_file.getline(buffer, 256);
			potential_file >> AB_p;
			potential_file.getline(buffer, 256);
			potential_file >> AB_q;
			potential_file.getline(buffer, 256);
			potential_file >> AB_r0;
			potential_file.getline(buffer, 256);
			potential_file >> BB_A;
			potential_file.getline(buffer, 256);
			potential_file >> BB_zeta;
			potential_file.getline(buffer, 256);
			potential_file >> BB_p;
			potential_file.getline(buffer, 256);
			potential_file >> BB_q;
			potential_file.getline(buffer, 256);
			potential_file >> BB_r0;
			potential_file.getline(buffer, 256);
			potential_file >> element_a;
			potential_file >> element_b;
		}
		else
		{
			potential_file >> element_a;
		}
	}

	potential_file.close();
}
