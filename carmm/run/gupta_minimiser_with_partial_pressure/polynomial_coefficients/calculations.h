#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include <math.h>
#include <stdio.h>

/**
@author Andrew Logsdail
*/
class Calculations{

public:
  Calculations();
  ~Calculations();
  void setNumbers(double AA, double pp, double qq, double xixi, double r0r0, double CsCs, double CeCe);
  double geta3() {return a3;}
  double geta4() {return a4;}
  double geta5() {return a5;}
  double getx3() {return x3;}
  double getx4() {return x4;}
  double getx5() {return x5;}
  void writeFile();
private:
  void calculate();
  double A, p, q, xi, r0, Cs, Ce;
  double AA_Cs, AA_Ce, AA_a3, AA_a4, AA_a5, AA_x3, AA_x4, AA_x5;
  double AB_Cs, AB_Ce, AB_a3, AB_a4, AB_a5, AB_x3, AB_x4, AB_x5;
  double a3, a4, a5, x3, x4, x5;
  int counter;
};

#endif
