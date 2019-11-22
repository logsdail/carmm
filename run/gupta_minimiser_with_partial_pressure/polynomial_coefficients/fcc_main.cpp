#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include "fcc_neighbours.h"

using namespace std;

int main(int argc, char *argv[])
{

  if(argc != 2)
  {
    //we don't have right number of arguments
    cout << "Catastrophic error: No r0 value provided" << endl;

    return EXIT_SUCCESS;
  }

  // Initial spacing so it looks pretty :)
  // cout << endl;
  // Now interpret parameters

  double r0 = 0;

  stringstream ss4(argv[1]);
  ss4 >> r0;

  // Addition to program so we use classes...
  // First off - cell neighbour lengths

  Fcc_neighbours fcc;
  fcc.setSides(r0);

  printf("===================\n");
  printf("NEIGHBOUR DISTANCES\n");
  printf("===================\n");

  printf ("Bond Length (r0) = %.16f \n", r0);
  printf ("First Neighbour = %.16f \n", fcc.getDistance(1.0));
  printf ("Second Neighbour = %.16f \n", fcc.getDistance(2.0));
  printf ("Third Neighbour = %.16f \n", fcc.getDistance(3.0));
  printf ("Fourth Neighbour = %.16f \n", fcc.getDistance(4.0));
  printf ("Fifth Neighbour = %.16f \n", fcc.getDistance(5.0));

}
