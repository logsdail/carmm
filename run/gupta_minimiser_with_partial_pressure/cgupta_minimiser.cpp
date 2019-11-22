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

CGupta_minimiser::CGupta_minimiser(string potential_filename, bool bCutoff, bool bVerbose)
{
	m_bCutoff = bCutoff; // Defines if the Cutoff is implemented
	m_bVerbose = bVerbose; // Defines if we are verbose in our output for atomic indices
	fxdx_counter = 0;

	// These are going to hold all the data from now on
        r0ij = new double[3];
        rhoij = new double[3];
        aaij = new double[3];
        z2ij = new double[3];
        qqij = new double[3];
        cc1ij = new double[3];
        cc2ij = new double[3];

	m_potentialFilename = potential_filename;
	this->ReadPotential();

	if (m_bCutoff){
                cutoff_start = new double[3];
		cutoff_end = new double[3];
		a3 = new double[3];
		a4 = new double[3];
		a5 = new double[3];
		x3 = new double[3];
		x4 = new double[3];
		x5 = new double[3];

        	this->ReadCutoff();
        }

        for (int i = 0; i < 3; i++)
	{
		z2ij[i] = z2ij[i]*z2ij[i];
		r0ij[i] = 1/r0ij[i];
		cc1ij[i] = -rhoij[i]*r0ij[i];
		cc2ij[i] = qqij[i]*r0ij[i];
		qqij[i] = 2*qqij[i];
	}
}

CGupta_minimiser::~CGupta_minimiser()
{
        delete [] r0ij;
        delete [] rhoij;
        delete [] aaij;
        delete [] z2ij;
        delete [] qqij;
        delete [] cc1ij;
        delete [] cc2ij;

        if (m_bCutoff){
        	delete [] cutoff_start;
	        delete [] cutoff_end;
        	delete [] a3;
	        delete [] a4;
        	delete [] a5;
	        delete [] x3;
        	delete [] x4;
	        delete [] x5;
	}
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
        static int nmax = n_atoms;
	double function_eval = 0.0;
        int ref = 0;

	// Increment our nice aesthetic counter
	fxdx_counter++;

	double v1i = 0.0;
	double v2i = 0.0;
	double v1ij = 0.0;
	double v2ij = 0.0;

	//shall place all of these on to the heap
	//so pointers away

	double *v1 = new double[nmax];
	double *v2 = new double[nmax];
	double *v2n = new double[nmax];

        double *dvdx = new double[nmax];
        double *dvdy = new double[nmax];
        double *dvdz = new double[nmax];

	int nmax_sqr_div = int((nmax*nmax)*0.5);

        double *dxij = new double[nmax_sqr_div];
        double *dyij = new double[nmax_sqr_div];
        double *dzij = new double[nmax_sqr_div];
        double *rij = new double[nmax_sqr_div];

        double *dv1dr = new double[nmax_sqr_div];
        double *dv2dr = new double[nmax_sqr_div];

        double *partial_energy = new double[nmax];
        double *partial_derivative = new double[nmax];
        double *partial_pressure = new double[nmax];

	double *volume = new double[nmax];

        int nn = 0;

	for(int i = 0; i < n_atoms; i++)
	{
		v1[i] = 0.0;
		v2[i] = 0.0;

                dvdx[i] = 0.0;
                dvdy[i] = 0.0;
                dvdz[i] = 0.0;

                // Partials
                partial_energy[i] = 0.0;
                partial_derivative[i] = 0.0;
                partial_pressure[i] = 0.0;

		// Volume assumes fcc
		// Take the determinant of the 3 x 3 matrix for the volume
		// Firstly, get the matrix values from bulk radius
		// As r = sqrt(x^2 + y^2 + z^2), this can be related to x = r/sqrt(2) for fcc

		if(m_Cluster.atom[i].atomType == element_a)
		{
			volume[i] = (1/r0ij[0])/sqrt(2);
		}
		else
		{
			volume[i] = (1/r0ij[2])/sqrt(2);
		}

		// Next up, take the determinant
		// By the book, this is a1b2c3 + a2b3c1 + a3b1c2 - c1b2a3 - cb3a1 - c3b1a2
		// But the diagonals are empty for fcc, so this cancels out to just a2b3c1 + a3b1c2
		// And as these sides are all the same, it becomes a cubic sum.

		volume[i] = 2.0 * (volume[i] * volume[i] * volume[i]);

                for(int j = i+1; j < n_atoms; j++)
		{
			nn++;
                        dxij[nn] = x[i] - x[j];
                        dyij[nn] = x[i+n_atoms] - x[j+n_atoms];
                        dzij[nn] = x[i+n_atoms*2] - x[j+n_atoms*2];
                        rij[nn] = sqrt(dxij[nn]*dxij[nn]+dyij[nn]*dyij[nn]+dzij[nn]*dzij[nn]);

                	// Need all this information for derivatives
	                dv1dr[nn] = 0.0;
        	        dv2dr[nn] = 0.0;
		}
	}

	nn = 0;

	for(int i = 0; i < n_atoms; i++) 
	{
		v1i = v1[i];
		v2i = v2[i];

		for(int j = i+1; j < n_atoms; j++)
		{
			nn++;
		
			// Reference atom type
			ref = 1;

			if(m_Cluster.atom[i].atomType == element_a && m_Cluster.atom[j].atomType == element_a)
			{
				ref = 0;
			} 
			else if(m_Cluster.atom[i].atomType == element_b && m_Cluster.atom[j].atomType == element_b)
			{
				ref = 2;
			}

			if (!m_bCutoff || rij[nn] < cutoff_start[ref]) 
			{
				v1ij = aaij[ref] * exp(-rhoij[ref] * ((rij[nn]*r0ij[ref]) - 1.0));
				v2ij = z2ij[ref] * exp(-qqij[ref] * ((rij[nn]*r0ij[ref]) - 1.0));

				dv1dr[nn] = 2 * v1ij * cc1ij[ref]; // Additional 2 for two-bodied component
				dv2dr[nn] = v2ij * cc2ij[ref];

			}
			else if (rij[nn] < cutoff_end[ref]) 
			{
	                        double factor1, factor2, factor3, factor4, factor5;

                                factor1 = rij[nn]-cutoff_end[ref];
                                factor2 = factor1*factor1;
                                factor3 = factor2*factor1;
                                factor4 = factor3*factor1;
                                factor5 = factor4*factor1;

                                v1ij = a5[ref]*factor5 + a4[ref]*factor4 + a3[ref]*factor3;
                                v2ij = x5[ref]*factor5 + x4[ref]*factor4 + x3[ref]*factor3;

				dv2dr[nn] = - ((5*x5[ref]*factor4 + 4*x4[ref]*factor3 + 3*x3[ref]*factor2)*v2ij);
				dv1dr[nn] = (5*a5[ref]*factor4 + 4*a4[ref]*factor3 + 3*a3[ref]*factor2);

				v2ij = v2ij*v2ij;
			} 
			else
			{
				
                                v1ij = 0;
                                v2ij = 0;

			}

	// End of current Andrew Logsdail additions....

			v1i += v1ij;
			v2i += v2ij;

			v1[j] += v1ij;
			v2[j] += v2ij;

		}

		v2n[i] = 1.0/sqrt(v2i);

		function_eval += (v1i - sqrt(v2i));
              
                // Save partial energy
                partial_energy[i] = (v1i -sqrt(v2i));
	}

	nn = 0;

	for(int i = 0; i < n_atoms; i++)
	{
		for(int j = i+1; j < n_atoms; j++)
		{
			nn++;

			// Note to self: This is two body, so we have 2*v1 + 2*v2 in the gradient
			// We just multiplied v1 by 2, but for v2 this depends on both 
			// v2n[i] and v2n[j], hence the sum is longer below.

                        dvdx[i] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dxij[nn] / rij[nn]);
                        dvdy[i] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dyij[nn] / rij[nn]);
                        dvdz[i] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dzij[nn] / rij[nn]);

			dvdx[j] -= (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dxij[nn] / rij[nn]);
			dvdy[j] -= (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dyij[nn] / rij[nn]);
			dvdz[j] -= (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j])) * (dzij[nn] / rij[nn]);

                        partial_derivative[i] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j]));
			partial_derivative[j] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j]));

                        partial_pressure[i] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j]))*rij[nn];
                        partial_pressure[j] += (dv1dr[nn] + dv2dr[nn] * (v2n[i] + v2n[j]))*rij[nn];
		}
	}

	//std::cout.precision(15);
	if (m_bVerbose)
	{
		std::cout << std::endl;
        	std::cout << "=========================================================================" << std::endl;
	        std::cout << " Fxdx Iteration = " << fxdx_counter << ", Energy = " << function_eval << " eV " << std::endl;
		std::cout << "=========================================================================" << std::endl;
	        std::cout << "     #     Type         E_{i}     dE_{i}/dr   (Numerical)         P_{i}" << std::endl;
        	std::cout << "=========================================================================" << std::endl;
	}

	for(int i = 0; i < n_atoms; i++)
	{
		// Partial pressures need scaling by a constant
		partial_pressure[i] *= (-1/(3*volume[i]));

		// Set gradients for L-BFGS
		g[i] = dvdx[i];
		g[i+n_atoms] = dvdy[i];
		g[i+2*n_atoms] = dvdz[i];

                // std::cout << g[i] << " " << g[i+n_atoms] << " " << g[i+2*n_atoms] << endl;

		if (m_bVerbose)
		{
			std::cout << setw(6) << i << setw(9) << m_Cluster.atom[i].atomType << setw(14) << partial_energy[i] << setw(14) << partial_derivative[i];
			std::cout << setw(14) << sqrt(g[i]*g[i]+g[i+n_atoms]*g[i+n_atoms]+g[i+2*n_atoms]*g[i+2*n_atoms]) << setw(14) << partial_pressure[i] << std::endl;
		}
	} 
        
	if (m_bVerbose)
	{
	        std::cout << "=========================================================================" << std::endl;
	}
	
	//ok we dont everything best delete the variables on the stack

	delete [] v1; 
	delete [] v2;
	delete [] v2n;
	delete [] dvdx;  
	delete [] dvdy;  
	delete [] dvdz;

        delete [] partial_energy;
        delete [] partial_derivative;
        delete [] partial_pressure;

	delete [] volume;

	delete [] dxij;
	delete [] dyij;
        delete [] dzij;
        delete [] rij;

        delete [] dv1dr;
        delete [] dv2dr;

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
                cutoff >> cutoff_start[0];
                cutoff.getline(buffer, 256);
                cutoff >> cutoff_end[0];
                cutoff.getline(buffer, 256);
                cutoff >> a3[0];
                cutoff.getline(buffer, 256);
                cutoff >> a4[0];
                cutoff.getline(buffer, 256);
                cutoff >> a5[0];
                cutoff.getline(buffer, 256);
                cutoff >> x3[0];
                cutoff.getline(buffer, 256);
                cutoff >> x4[0];
                cutoff.getline(buffer, 256);
                cutoff >> x5[0];
                cutoff.getline(buffer, 256);

		if ( num_els_cutoff == 2 )
		{
                        cutoff >> cutoff_start[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> cutoff_end[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> a3[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> a4[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> a5[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> x3[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> x4[1];
                        cutoff.getline(buffer, 256);
                        cutoff >> x5[1];
                        cutoff.getline(buffer, 256);

			cutoff >> cutoff_start[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> cutoff_end[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> a3[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> a4[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> a5[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> x3[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> x4[2];
                	cutoff.getline(buffer, 256);
                	cutoff >> x5[2];
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
		potential_file >> aaij[0];
		potential_file.getline(buffer, 256);
		potential_file >> z2ij[0];
		potential_file.getline(buffer, 256);
		potential_file >> rhoij[0];
		potential_file.getline(buffer, 256);
		potential_file >> qqij[0];
		potential_file.getline(buffer, 256);
		potential_file >> r0ij[0];
		potential_file.getline(buffer, 256);

		if ( num_els == 2 )
		{
			potential_file >> aaij[1];
			potential_file.getline(buffer, 256);
			potential_file >> z2ij[1];
			potential_file.getline(buffer, 256);
			potential_file >> rhoij[1];
			potential_file.getline(buffer, 256);
			potential_file >> qqij[1];
			potential_file.getline(buffer, 256);
			potential_file >> r0ij[1];
			potential_file.getline(buffer, 256);
			potential_file >> aaij[2];
			potential_file.getline(buffer, 256);
			potential_file >> z2ij[2];
			potential_file.getline(buffer, 256);
			potential_file >> rhoij[2];
			potential_file.getline(buffer, 256);
			potential_file >> qqij[2];
			potential_file.getline(buffer, 256);
			potential_file >> r0ij[2];
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
