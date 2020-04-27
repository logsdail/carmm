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
#include "cminimiser.h"
#include "Struct.h"
#include "cError.h"

using namespace std;
//Blank constructor intialising with methods required
CMinimiser::CMinimiser()
{

};

//bad why impose problem specific needs on minimiser?
//Constructor uses existing vector of atoms
CMinimiser::~CMinimiser()
{
}

//Setups and runs the lbfgs minimiser algorithm uses the pure abstract function fx_dx which is implemented in dervied classes ie CGupta_Minimiser
//Changed to use scitbx::lbfgs::minimizer()
void CMinimiser::minimise()
{
	//nmax = number of variables plus one
	int nmax = m_inumVariables + 1; //this is how chris sets his array(take max mem) this must be 3 times bigger than the number of atoms in a system
	int mmax = 20;

	//n = number of variables
	int n = m_inumVariables; //was 3 * n_atoms
	int m = 5;		//suitable value
	//again deprecated or at least has anaologues in minimizer constructer
	double factr = 1.0e5;
	double pgtol = 1.0e-7;

	//we will try first with just initalising the minimiser using most defaults
	//place on the heap
	m_pMinimizer = new scitbx::lbfgs::minimizer<double>(n, m, 20);

	//also we want a convergence check we use the one that comes with this library
	scitbx::lbfgs::traditional_convergence_test<double> is_Converged(n);

	//minimizer is all ready to go so we shall check our settings
	this->PrintMinimiserOptions();

	//create to standard vectors with all to place all our variables
	std::vector<double> x(3*nmax); //nmax = variable nums, this is to much it should just be all the variables
	std::vector<double> g(3*nmax); //same i will change this after I have it working like before.

	//function evaluation.
	double f = 0.0;

	//pointers to the first element of variable and gradient vectors
	double* ptrX = &(*(x.begin()));
	double* ptrG = &(*(g.begin()));

	for(int i = 0; i < m_inumVariables; i++)
	{
		x[i] = m_vVariables[i]; //zero out and assign variables however they both vectors now so this is a bit ott
		g[i] = 0.0;
	}

	for(;;)//an infinite for loop so must have break case
	{
		//call the fx_dx implemented by the concrete class
		f = this->fx_dx(ptrX,ptrG);

		//if minimizer returns true it requires fx_dx calc so restart the loop
		if(m_pMinimizer->run(ptrX, f, ptrG))
		{
			std::cout << "fxdx requested" << std::endl;
			continue;
		}

		//if it returns false we should check for convergence
		if(is_Converged(ptrX,ptrG))
		{
			std::cout << "Convergence obtained" << std::endl;
			std::cout << "gnorm" << m_pMinimizer->euclidean_norm(&(*(g.begin()))) << std::endl;
			std::cout << "nfunc: " << m_pMinimizer->nfun() << std::endl;
			std::cout << "iter: " << m_pMinimizer->iter() << std::endl;
			break;
		}
		//if we are doing excessive function evaluations
		if(m_pMinimizer->nfun() > 2000)
		{
			std::cout << "Warning large number of fucntion evaluations occuring exiting minimiser\n";
			break;
		}

		m_pMinimizer->run(ptrX, f, ptrG);
	}	

	m_dEnergy = f;
	for(int i = 0; i < m_inumVariables; i++)
	{
		m_vVariables[i] = x[i];
	}
	//ok whats going on here? 
		//release our dynamic objects

	//fine to do this
	delete m_pMinimizer;
		//these are both pointers to the first value of a vector
	//vector will handle itself so I can safely derefence them
	ptrG = NULL;
	ptrX = NULL;
}

void CMinimiser::CalculateEnergy()
{
	int nmax = m_inumVariables + 1; //this is how chris sets his array(take max mem) this must be 3 times bigger than the number of atoms in a system

	//create to standard vectors with all to place all our variables
	std::vector<double> x(3*nmax); //nmax = variable nums, this is to much it should just be all the variables
	std::vector<double> g(3*nmax); //same i will change this after I have it working like before.

	//function evaluation.
	double f = 0.0;

	//pointers to the first element of variable and gradient vectors
	double* ptrX = &(*(x.begin()));
	double* ptrG = &(*(g.begin()));

	for(int i = 0; i < m_inumVariables; i++)
	{
		x[i] = m_vVariables[i]; //zero out and assign variables however they both vectors now so this is a bit ott
		g[i] = 0.0;
	}

	m_dEnergy = this->fx_dx(ptrX, ptrG);
}


void CMinimiser::PrintMinimiserOptions()
{
	if(!m_pMinimizer)
	{
		throw cError("No minimiser instantiated!");
	}
	std::cout << "======= Minimiser Settings =======\n";
	std::cout << "n: " << m_pMinimizer->n() << "\n";
	std::cout << "m: " << m_pMinimizer->m() << "\n";
	std::cout << "xtol: " << m_pMinimizer->xtol() << "\n";
	std::cout << "gtol: " << m_pMinimizer->gtol() << "\n";
	std::cout << "stpmin: " << m_pMinimizer->stpmin() << "\n";
	std::cout << "stpmax: " << m_pMinimizer->stpmax() << "\n";
}
