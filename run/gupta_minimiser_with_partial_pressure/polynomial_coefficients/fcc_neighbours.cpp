#include <math.h>
#include "fcc_neighbours.h"

double n1;

Fcc_neighbours::Fcc_neighbours()
{
  n1 = 0;
}

Fcc_neighbours::~Fcc_neighbours() {}

void Fcc_neighbours::setSides(const double input)
{
  n1 = input;
}

double Fcc_neighbours::getDistance(const double input)
{
  return (n1 * sqrt (input));
}
