#ifndef FCC_NEIGHBOURS_H
#define FCC_NEIGHBOURS_H
#include <math.h>

/**
@author Andrew Logsdail
*/
class Fcc_neighbours{

public:
  Fcc_neighbours();
  ~Fcc_neighbours();
  void setSides(double input);
  double getDistance(double input);
private:
  double n1;

};

#endif
