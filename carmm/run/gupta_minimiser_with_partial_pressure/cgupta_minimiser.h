/***************************************************************************
 *   Copyright (C) 2005 by Ben Curley                                      *
 *   curley@tc.bham.ac.uk                                                  *
 *	

	Complete and utter crap, I gonna hack a solution but god what was I doing?
	This class should not care or know of the problem it solves its only concern
	is the variables and gradient function, but here i preforming conversions
	file reading and other gubbings for no good reason. Same for cminimiser
	which this class inherits. I mean come on inherting from an abstract base 
    class a file reading method! 
 ***************************************************************************/
#ifndef CGUPTA_MINIMISER_H
#define CGUPTA_MINIMISER_H

#include "cminimiser.h"
#include "Struct.h"
#include <string>
/**
@author Ben Curley
*/
class CGupta_minimiser : public CMinimiser
{
public:

	CGupta_minimiser(string potential_filename, bool bCutoff, bool bVerbose);
	virtual ~CGupta_minimiser();

	virtual double fx_dx(const double* x, double* g); //Derivation of gupta pot called by cminimiser

	void SetClusterData(clusterData cluster);
	clusterData GetClusterData();
	clusterData GetMinimisedData(); //Gaurentees minimisation
	double GetEnergy();
private:
	void ReadPotential(); //not best place but it is ok
	void ReadCutoff(); // I don't well enough if this is the best place or not but it works :)
	clusterData m_Cluster;
	std::string m_potentialFilename;

	int num_els, num_els_cutoff;
	bool m_bCutoff; // Implies cutoff will be used. Cutoff file must be present
	bool m_bVerbose; // Defines verbosity of output

	string element_a;
	string element_b;

        double *r0ij;
        double *rhoij;
        double *aaij;
        double *z2ij;
        double *qqij;
        double *cc1ij;
        double *cc2ij;

	double *cutoff_start;
	double *cutoff_end;
	double *a3;
	double *a4;
	double *a5;
	double *x3;
	double *x4;
	double *x5;

	int fxdx_counter;
};

#endif
