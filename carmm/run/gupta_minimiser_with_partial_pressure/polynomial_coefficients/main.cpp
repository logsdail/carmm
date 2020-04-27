#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fcc_neighbours.h"
#include "calculations.h"

using namespace std;

int main(int argc, char *argv[])
{

  if(argc != 4)
  {
    //we don't have right number of arguments
    cout << "Invalid number of arguments" << endl;
    cout << "You need, in this order: Potential file, Cutoff start and Cutoff end (in relation to neighbours)" << endl;
    cout << "i.e. ./polynomial gupta_Au.in 2 3 will give you a cutoff between second and third neighbours. Input does allow doubles for cutoff positions." << endl;
    return EXIT_SUCCESS;
  }

  // Initial spacing so it looks pretty :)
  // cout << endl;
  // Now interpret parameters

  double Cs, Ce = 0.0; 
  int num_els = 0;
  double AA_A, AA_p, AA_q, AA_xi, AA_r0, AA_Cs, AA_Ce = 0;
  double AB_A, AB_p, AB_q, AB_xi, AB_r0, AB_Cs, AB_Ce = 0;
  double BB_A, BB_p, BB_q, BB_xi, BB_r0, BB_Cs, BB_Ce = 0;

  string parameters = argv[1];
  fstream potential_file;
  char buffer[256];

  potential_file.open(parameters.c_str(), ios::in);
  if(!potential_file.is_open())
    return 1;
                
  potential_file >> num_els;
  potential_file >> AA_A;
  potential_file.getline(buffer, 256);
  potential_file >> AA_xi;
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
      potential_file >> AB_xi;
      potential_file.getline(buffer, 256);
      potential_file >> AB_p;
      potential_file.getline(buffer, 256);
      potential_file >> AB_q;
      potential_file.getline(buffer, 256);
      potential_file >> AB_r0;
      potential_file.getline(buffer, 256);
      potential_file >> BB_A;
      potential_file.getline(buffer, 256);
      potential_file >> BB_xi;
      potential_file.getline(buffer, 256);
      potential_file >> BB_p;
      potential_file.getline(buffer, 256);
      potential_file >> BB_q;
      potential_file.getline(buffer, 256);
      potential_file >> BB_r0;
      potential_file.getline(buffer, 256);
    }
  potential_file.close(); 

  //stringstream ss5(argv[2]);
  //ss5 >> Cs;

  //stringstream ss6(argv[3]);
  //ss6 >> Ce;

  Cs = atof(argv[2]);
  Ce = atof(argv[3]);

/////////////////////////////////////////////////////////////////////////////////////

  // First off - cell neighbour lengths

  Fcc_neighbours fcc;
  fcc.setSides(AA_r0);
  AA_Cs = fcc.getDistance(Cs);
  AA_Ce = fcc.getDistance(Ce);
  if ( num_els == 2 )
  {
    fcc.setSides(AB_r0);
    AB_Cs = fcc.getDistance(Cs);
    AB_Ce = fcc.getDistance(Ce);
    fcc.setSides(BB_r0);
    BB_Cs = fcc.getDistance(Cs);
    BB_Ce = fcc.getDistance(Ce);    
  }

//////////////////////////////////////////////////////////////////////////////////////

  // Now calculations

  Calculations calc;
  calc.setNumbers(AA_A,AA_p,AA_q,AA_xi,AA_r0,AA_Cs,AA_Ce);

  if ( num_els == 2 )
    {
      calc.setNumbers(AB_A,AB_p,AB_q,AB_xi,AB_r0,AB_Cs,AB_Ce);
      calc.setNumbers(BB_A,BB_p,BB_q,BB_xi,BB_r0,BB_Cs,BB_Ce);
    }

//////////////////////////////////////////////////////////////////////////////////////

  // And finish by writing to file

  calc.writeFile();

return EXIT_SUCCESS;
}

