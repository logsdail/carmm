/*
	CMinimiser: Base Abstract class which handles the LBFGS wrapper

*/
#ifndef CMINIMISER_H
#define CMINIMISER_H

#include "lbfgs.h"
#include "Struct.h"

class CMinimiser
{
public:
	//constructor 
	CMinimiser();
	virtual ~CMinimiser();
	virtual void minimise(); //Call this to preform minimisation
	virtual double fx_dx(const double* x, double* g) = 0; //Must be implemented in concrete classes

	void CalculateEnergy();									   //Performs fxdx on current var.
	void SetNumVars(int iNumVars){m_inumVariables = iNumVars;} //Set size of problem
	int GetNumVars(){return m_inumVariables;}				   //Get size of problem
	void SetVarVector(vector<double> vec){m_vVariables = vec;} //Set Problem set
	vector<double> GetVarVector(){return m_vVariables;}		   //Get Problem set
	void PrintMinimiserOptions();
protected:
	int m_inumVariables; //number of variables in the problem. i.e for 1 atom x y z = 3 variables
	double m_dEnergy; //minimisation energy based on potential used.
	vector<double> m_vVariables; //vector holding variables
	scitbx::lbfgs::minimizer<double>* m_pMinimizer;
};


	/*
class CMinimiser{
public:
	CMinimiser();
	CMinimiser(vector<atom> vec);
	CMinimiser(string structure_filename);
	virtual ~CMinimiser();
	virtual void minimise(); 
	virtual double fx_dx(const double x[], double g[]) = 0; //Pure virtual method required to compile this base class in now abstract
	void set_n_atoms(int x){n_atoms = x;}
	void set_atom_vec(vector<atom> vec){atom_vec = vec;}
	int get_n_atoms(){return n_atoms;}
	double get_energy(){return energy;}
	vector<atom> get_atom_vec(){return atom_vec;}
protected:
	int n_atoms;
	vector<atom> atom_vec;
	double energy;
};
*/
#endif
