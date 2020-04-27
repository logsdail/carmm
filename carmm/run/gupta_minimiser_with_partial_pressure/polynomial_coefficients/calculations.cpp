#include <stdio.h>
#include <math.h>
#include "calculations.h"

using namespace std;

double A, p, q, xi, r0, Cs, Ce;
double AA_Cs, AA_Ce, AA_a3, AA_a4, AA_a5, AA_x3, AA_x4, AA_x5;
double AB_Cs, AB_Ce, AB_a3, AB_a4, AB_a5, AB_x3, AB_x4, AB_x5;
double a3, a4, a5, x3, x4, x5;
int counter;

Calculations::Calculations()
{
  counter = 0;
  A = 0;
  p = 0;
  q = 0;
  xi = 0;
  r0 = 0;
  Cs = 0;
  Ce = 0;
  a3 = 0;
  a4 = 0;
  a5 = 0;
  x3 = 0;
  x4 = 0;
  x5 = 0;
}

Calculations::~Calculations() {}

void Calculations::setNumbers(double AA, double pp, double qq, double xixi, double r0r0, double CsCs, double CeCe)
{
  if ( counter == 1 )
  {
    AA_Cs = Cs;
    AA_Ce = Ce;
    AA_a3 = a3;
    AA_a4 = a4;
    AA_a5 = a5;
    AA_x3 = x3;
    AA_x4 = x4;
    AA_x5 = x5;
  } 
  else if ( counter == 2 )
  {
    AB_Cs = Cs;
    AB_Ce = Ce;
    AB_a3 = a3;
    AB_a4 = a4;
    AB_a5 = a5;
    AB_x3 = x3;
    AB_x4 = x4;
    AB_x5 = x5;
  }

  A = AA;
  p = pp;
  q = qq;
  xi = xixi;
  r0 = r0r0;
  Cs = CsCs;
  Ce = CeCe;

  a3 = 0;
  a4 = 0;
  a5 = 0;
  x3 = 0;
  x4 = 0;
  x5 = 0;

  counter++;

//  cout << a << endl;
//  cout << p << endl;
//  cout << q << endl;
//  cout << xi << endl;
//  cout << r0 << endl;

  calculate();
}

void Calculations::calculate()
{

  double rDivr0Minus1 = (Cs/r0) - 1;

  double rMinusCe, rMinusCe2, rMinusCe3, rMinusCe4, rMinusCe5;
  rMinusCe = Cs - Ce;
  rMinusCe2 = rMinusCe * rMinusCe;
  rMinusCe3 = rMinusCe * rMinusCe2;
  rMinusCe4 = rMinusCe * rMinusCe3;
  rMinusCe5 = rMinusCe * rMinusCe4;

  double FiveTimesrMinusCe4, FourTimesrMinusCe3, ThreeTimesrMinusCe2;
  FiveTimesrMinusCe4 = 5 * rMinusCe4;
  FourTimesrMinusCe3 = 4 * rMinusCe3;
  ThreeTimesrMinusCe2 = 3 * rMinusCe2;

  double TwentyTimesrMinusCe3, TwelveTimesrMinusCe2, SixTimesrMinusCe;
  TwentyTimesrMinusCe3 = 20 * rMinusCe3;
  TwelveTimesrMinusCe2 = 12 * rMinusCe2;
  SixTimesrMinusCe = 6 * rMinusCe;

  // Sort out variables for Repulsive term
  double MinuspTimesrDivExpression, AexpMinuspExpression = 0;
  MinuspTimesrDivExpression = -p * rDivr0Minus1; //-p((r/r0)-1)
  AexpMinuspExpression = A * exp(MinuspTimesrDivExpression); //A*exp[-p((r/r0)-1)]

  // Firstly, A*exp[-p((r/r0)-1)] = a5(r-Ce)^5 + a4(r-Ce)^4 + a3(r-Ce)^3
  //cout << "Start:  " << AexpMinuspExpression << " = " << rMinusCe5 << "(a5) + " << rMinusCe4 << "(a4) + " << rMinusCe3 << "(a3)" <<  endl;
  // Next, (-p/r0)A*exp[-p((r/r0)-1)] = 5*a5(r-Ce)^4 + 4*a4(r-Ce)^3 + 3*a3(r-Ce)^2

  double FirstDerivativeAexpMinuspExpression = 0;
  FirstDerivativeAexpMinuspExpression = ( -p / r0) * AexpMinuspExpression; //(-p/r0)*A*exp[-p((r/r0)-1)]

  //cout << "First Derivative:  " << FirstDerivativeAexpMinuspExpression << " = " << FiveTimesrMinusCe4 << "(a5) + " << FourTimesrMinusCe3 << "(a4) + " << ThreeTimesrMinusCe2 << "(a3)" << endl;
  // Finally for repulsive part, Second Derivatives - (-p/r0)^2*A*exp[-p((r/r0)-1)] = 20*a5(r-Ce)^3 + 12*a4(r-Ce)^2 + 6*A3(r-Ce)

  double SecondDerivativeAexpMinuspExpression = 0;
  SecondDerivativeAexpMinuspExpression = ( -p / r0) * FirstDerivativeAexpMinuspExpression; // (-p/r0)^2*A*exp[-p((r/r0)-1)]

  // cout << "Second Derivative:  " << SecondDerivativeAexpMinuspExpression << " = " << TwentyTimesrMinusCe3 << "(a5) + " << TwelveTimesrMinusCe2 << "(a4) + " << SixTimesrMinusCe << "(a3)" << endl;

  // Now lets work out a5, a4 and a3. Need to narrow this down to a bi linear equation sotime to equal out.
  // Will call this first two variable equation number 1, and the second number 2, and name the variables to match

  // Equation 1
  // First Derivative vs. Original

  double Eq1Factor, Eq1OriginalRepulsive, Eq1Originala5, Eq1Originala4, Eq1Originala3 = 0;
  Eq1Factor = - (FiveTimesrMinusCe4 / rMinusCe5);
  Eq1OriginalRepulsive = AexpMinuspExpression * Eq1Factor;
  Eq1Originala5 = rMinusCe5 * Eq1Factor;
  Eq1Originala4 = rMinusCe4 * Eq1Factor;
  Eq1Originala3 = rMinusCe3 * Eq1Factor;

  double Eq1Repulsive, Eq1a5, Eq1a4, Eq1a3 = 0;
  Eq1Repulsive = Eq1OriginalRepulsive + FirstDerivativeAexpMinuspExpression;
  Eq1a5 = Eq1Originala5 + FiveTimesrMinusCe4;
  Eq1a4 = Eq1Originala4 + FourTimesrMinusCe3;
  Eq1a3 = Eq1Originala3 + ThreeTimesrMinusCe2;

  // Equation 2
  // FirstDerivative vs. Second Derivative
  double Eq2Factor, Eq2FirstRepulsive, Eq2Firsta5, Eq2Firsta4, Eq2Firsta3 = 0;
  Eq2Factor = - (TwentyTimesrMinusCe3 / FiveTimesrMinusCe4);
  Eq2FirstRepulsive = FirstDerivativeAexpMinuspExpression * Eq2Factor;
  Eq2Firsta5 = FiveTimesrMinusCe4 * Eq2Factor;
  Eq2Firsta4 = FourTimesrMinusCe3 * Eq2Factor;
  Eq2Firsta3 = ThreeTimesrMinusCe2 * Eq2Factor;

  double Eq2Repulsive, Eq2a5, Eq2a4, Eq2a3 = 0;
  Eq2Repulsive = Eq2FirstRepulsive + SecondDerivativeAexpMinuspExpression;
  Eq2a5 = Eq2Firsta5 + TwentyTimesrMinusCe3;
  Eq2a4 = Eq2Firsta4 + TwelveTimesrMinusCe2;
  Eq2a3 = Eq2Firsta3 + SixTimesrMinusCe;

  // So bilinear equations created - let's put these on the screen so we can check them
  // cout << endl;
  // cout << "Now to linear equations - here is what we have calculated:" << endl;
  // cout << "Equation 1: " << Eq1Repulsive << " = " << Eq1a5 << "(a5) + " << Eq1a4 << "(a4) + " << Eq1a3 << "(a3)" << endl;
  // cout << "Equation 2: " << Eq2Repulsive << " = " << Eq2a5 << "(a5) + " << Eq2a4 << "(a4) + " << Eq2a3 << "(a3)" << endl;
  // cout << "We can negate the a5 factor as it is >> 0" << endl;
  // cout << endl;

  // Now to solve the linear equations we have
  // Let's start by working out a3
  // Equation 3!

  double Eq3Factor, Eq3OriginalRepulsive, Eq3Originala4, Eq3Originala3 = 0;
  Eq3Factor = - (Eq2a4 / Eq1a4);
  Eq3OriginalRepulsive = Eq1Repulsive * Eq3Factor;
  Eq3Originala4 = Eq1a4 * Eq3Factor;
  Eq3Originala3 = Eq1a3 * Eq3Factor;

  // cout << "So now we have:" << endl;
  // cout << "Equation 3: " << Eq3OriginalRepulsive << " = " << Eq3Originala4 << "(a4) + " << Eq3Originala3 << "(a3)" << endl;
  // cout << "Equation 2: " << Eq2Repulsive << " = " << Eq2a4 << "(a4) + " << Eq2a3 << "(a3)" << endl;
  // cout << endl;

  double solvea3Repulsive, solvea3;
  solvea3Repulsive = Eq3OriginalRepulsive + Eq2Repulsive;
  solvea3 = Eq3Originala3 + Eq2a3;

  // cout << "Add these together and we get..." << endl;
  // cout << solvea3Repulsive << " = " << Eq3Originala4 + Eq2a4 << "(a4) +" << solvea3 << "(a3)" << endl;
  // cout << "which means..." << endl;

  a3 = solvea3Repulsive / solvea3;
  // printf ("a3 = %.16f \n", a3); // cout<<endl;

  a4 = (Eq2Repulsive - (Eq2a3 * a3)) / Eq2a4;
  // printf ("a4 = %.16f \n", a4); // cout<<endl;

  a5 = ((AexpMinuspExpression - (rMinusCe3 * a3) - (rMinusCe4 * a4)) / rMinusCe5);
  // printf ("a5 = %.16f \n", a5); 

  // Now for attractive terms!

  // cout << endl;
  // cout << "===============" << endl;
  // cout << "ATTRACTIVE TERM" << endl;
  // cout << "===============" << endl;

  // Sort out variables for attractive term
  double MinusqTimesrDivExpression, XiexpMinusqExpression = 0;
  MinusqTimesrDivExpression = -q * rDivr0Minus1; //-q((r/r0)-1)
  XiexpMinusqExpression = xi * exp(MinusqTimesrDivExpression); //xi*exp[-q((r/r0)-1)]

  // Firstly, xi*exp[-q((r/r0)-1)] = x5(r-Ce)^5 + x4(r-Ce)^4 + x3(r-Ce)^3
  // cout << "Start:  " << XiexpMinusqExpression << " = " << rMinusCe5 << "(x5) + " << rMinusCe4 << "(x4) + " << rMinusCe3 << "(x3)" <<  endl;

  // Next, (-q/r0)xi*exp[-q((r/r0)-1)] = 5*x5(r-Ce)^4 + 4*x4(r-Ce)^3 + 3*x3(r-Ce)^2
  double FirstDerivativeXiexpMinusqExpression = 0;
  FirstDerivativeXiexpMinusqExpression = ( -q / r0) * XiexpMinusqExpression; //(-q/r0)*xi*exp[-q((r/r0)-1)]

  // cout << "First Derivative:  " << FirstDerivativeXiexpMinusqExpression << " = " << FiveTimesrMinusCe4 << "(x5) + " << FourTimesrMinusCe3 << "(x4) + " << ThreeTimesrMinusCe2 << "(x3)" << endl;
  // Finally for repulsive part, Second Derivatives - (-q/r0)^2*xi*exp[-q((r/r0)-1)] = 20*x5(r-Ce)^3 + 12*x4(r-Ce)^2 + 6*x3(r-Ce)

  double SecondDerivativeXiexpMinusqExpression = 0;
  SecondDerivativeXiexpMinusqExpression = ( -q / r0) * FirstDerivativeXiexpMinusqExpression; // (-q/r0)^2*xi*exp[-q((r/r0)-1)]

  // cout << "Second Derivative:  " << SecondDerivativeXiexpMinusqExpression << " = " << TwentyTimesrMinusCe3 << "(x5) + " << TwelveTimesrMinusCe2 << "(x4) + " << SixTimesrMinusCe << "(x3)" << endl;
  // Now lets work out x5, x4 and x3. Need to narrow this down to a bi linear equation sotime to equal out.
  // Will call this first two variable equation number 4, and the second number 5, and name the variables to match

  // Equation 4
  // First Derivative vs. Original - THESE VARIABLES HAVE SAME VALUES AS FOR EQUATION 1

  double Eq4Factor, Eq4OriginalAttractive, Eq4Originalx5, Eq4Originalx4, Eq4Originalx3 = 0;
  Eq4Factor = - (FiveTimesrMinusCe4 / rMinusCe5);
  Eq4OriginalAttractive = XiexpMinusqExpression * Eq4Factor;
  Eq4Originalx5 = rMinusCe5 * Eq4Factor;
  Eq4Originalx4 = rMinusCe4 * Eq4Factor;
  Eq4Originalx3 = rMinusCe3 * Eq4Factor;

  double Eq4Attractive, Eq4x5, Eq4x4, Eq4x3 = 0;
  Eq4Attractive = Eq4OriginalAttractive + FirstDerivativeXiexpMinusqExpression;
  Eq4x5 = Eq4Originalx5 + FiveTimesrMinusCe4;
  Eq4x4 = Eq4Originalx4 + FourTimesrMinusCe3;
  Eq4x3 = Eq4Originalx3 + ThreeTimesrMinusCe2;

  // Equation 5
  // FirstDerivative vs. Second Derivative
  double Eq5Factor, Eq5FirstAttractive, Eq5Firstx5, Eq5Firstx4, Eq5Firstx3 = 0;
  Eq5Factor = - (TwentyTimesrMinusCe3 / FiveTimesrMinusCe4);
  Eq5FirstAttractive = FirstDerivativeXiexpMinusqExpression * Eq5Factor;
  Eq5Firstx5 = FiveTimesrMinusCe4 * Eq5Factor;
  Eq5Firstx4 = FourTimesrMinusCe3 * Eq5Factor;
  Eq5Firstx3 = ThreeTimesrMinusCe2 * Eq5Factor;

  double Eq5Attractive, Eq5x5, Eq5x4, Eq5x3 = 0;
  Eq5Attractive = Eq5FirstAttractive + SecondDerivativeXiexpMinusqExpression;
  Eq5x5 = Eq5Firstx5 + TwentyTimesrMinusCe3;
  Eq5x4 = Eq5Firstx4 + TwelveTimesrMinusCe2;
  Eq5x3 = Eq5Firstx3 + SixTimesrMinusCe;

  // So bilinear equations created - let's put these on the screen so we can check them
  // cout << endl;
  // cout << "Now to linear equations - here is what we have calculated:" << endl;
  // cout << "Equation 4: " << Eq4Attractive << " = " << Eq4x5 << "(x5) + " << Eq4x4 << "(x4) + " << Eq4x3 << "(x3)" << endl;
  // cout << "Equation 5: " << Eq5Attractive << " = " << Eq5x5 << "(x5) + " << Eq5x4 << "(x4) + " << Eq5x3 << "(x3)" << endl;
  // cout << "We can negate the x5 factor as it is >> 0" << endl;
  // cout << endl;

  // Now to solve the linear equations we have
  // Let's start by working out x3
  // Equation 6!

  double Eq6Factor, Eq6OriginalAttractive, Eq6Originalx4, Eq6Originalx3 = 0;
  Eq6Factor = - (Eq5x4 / Eq4x4);
  Eq6OriginalAttractive = Eq4Attractive * Eq6Factor;
  Eq6Originalx4 = Eq4x4 * Eq6Factor;
  Eq6Originalx3 = Eq4x3 * Eq6Factor;

  // cout << "So now we have:" << endl;
  // cout << "Equation 6: " << Eq6OriginalAttractive << " = " << Eq6Originalx4 << "(x4) + " << Eq6Originalx3 << "(x3)" << endl;
  // cout << "Equation 5: " << Eq5Attractive << " = " << Eq5x4 << "(x4) + " << Eq5x3 << "(x3)" << endl;
  // cout << endl;

  double solvex3Attractive, solvex3;
  solvex3Attractive = Eq6OriginalAttractive + Eq5Attractive;
  solvex3 = Eq6Originalx3 + Eq5x3;

  // cout << "Add these together and we get..." << endl;
  // cout << solvex3Attractive << " = " << Eq6Originalx4 + Eq5x4 << "(x4) +" << solvex3 << "(x3)" << endl;
  // cout << "which means..." << endl;

  x3 = solvex3Attractive / solvex3;
  // printf ("x3 = %.16f \n", x3);

  x4 = (Eq5Attractive - (Eq5x3 * x3)) / Eq5x4;
  // printf ("x4 = %.16f \n", x4);

  x5 = ((XiexpMinusqExpression - (rMinusCe3 * x3) - (rMinusCe4 * x4)) / rMinusCe5);
  // printf ("x5 = %.16f \n", x5);
  // cout << endl;
}

void Calculations::writeFile()
{
  FILE *outData;
  bool fIsOpen = ( outData = fopen("cutoff_parameters", "w") );

  if ( fIsOpen )
  {
    if (counter == 1)
    {
      fprintf(outData, "1\n");
    } 
    else
    {
      fprintf(outData, "2\n");
      fprintf(outData, "%.16f     AA cutoff_start \n", AA_Cs);
      fprintf(outData, "%.16f     AA cutoff_end   \n", AA_Ce);
      fprintf(outData, "%.16f     AA a3           \n", AA_a3);
      fprintf(outData, "%.16f     AA a4           \n", AA_a4);
      fprintf(outData, "%.16f     AA a5           \n", AA_a5);
      fprintf(outData, "%.16f     AA x3           \n", AA_x3);
      fprintf(outData, "%.16f     AA x4           \n", AA_x4);
      fprintf(outData, "%.16f     AA x5           \n", AA_x5);
      fprintf(outData, "%.16f     AB cutoff_start \n", AB_Cs);
      fprintf(outData, "%.16f     AB cutoff_end   \n", AB_Ce);
      fprintf(outData, "%.16f     AB a3           \n", AB_a3);
      fprintf(outData, "%.16f     AB a4           \n", AB_a4);
      fprintf(outData, "%.16f     AB a5           \n", AB_a5);
      fprintf(outData, "%.16f     AB x3           \n", AB_x3);
      fprintf(outData, "%.16f     AB x4           \n", AB_x4);
      fprintf(outData, "%.16f     AB x5           \n", AB_x5);
    }

    fprintf(outData, "%.16f     cutoff_start \n", Cs);
    fprintf(outData, "%.16f     cutoff_end   \n", Ce);
    fprintf(outData, "%.16f     a3           \n", a3);
    fprintf(outData, "%.16f     a4           \n", a4);
    fprintf(outData, "%.16f     a5           \n", a5);
    fprintf(outData, "%.16f     x3           \n", x3);
    fprintf(outData, "%.16f     x4           \n", x4);
    fprintf(outData, "%.16f     x5           \n", x5);
    fclose(outData);
  }
  else
    printf("Catastrophic error opening cutoff_parameters for writing\n");

}

