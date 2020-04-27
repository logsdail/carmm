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

	CGupta_minimiser(string potential_filename, bool bCutoff);
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

	string element_a;
	string element_b;

	double AA_r0;
	double AA_p;
	double AA_A;
	double AA_q;
	double AA_zeta;

	double AB_r0;
	double AB_p;
	double AB_A;
	double AB_q;
	double AB_zeta;

	double BB_r0;
	double BB_p;
	double BB_A;
	double BB_q;
	double BB_zeta;

	double AA_cutoff_start;
	double AA_cutoff_end;
	double AA_a3;
	double AA_a4;
	double AA_a5;
	double AA_x3;
	double AA_x4;
	double AA_x5;

        double BB_cutoff_start;
        double BB_cutoff_end;
        double BB_a3;
        double BB_a4;
        double BB_a5;
        double BB_x3;
        double BB_x4;
        double BB_x5;

        double AB_cutoff_start;
        double AB_cutoff_end;
        double AB_a3;
        double AB_a4;
        double AB_a5;
        double AB_x3;
        double AB_x4;
        double AB_x5;

};

#endif
