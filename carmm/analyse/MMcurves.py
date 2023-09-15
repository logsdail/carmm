# A few functions and scripts for generating and analysing interatomic separation curves for simple systems.
# Authored by Ahment Sayin as part of a Nuffield project placement during August 2023.

# Buckingham Potential
def buckingham_potential(r, A, rho, C):
  # r float: Interionic separation
  # A float: The A parameter for a Buckingham potential
  # rho float: The damping parameter on the separation
  # C float: Defines the strength of the dispersive forces on the ions
  import numpy as np
  
  buck = A * np.exp(-r / rho) - (C / r ** 6)
  
  return buck

# Lennard-Jones Potential
def lj_potential(r, e, E0, E, S):
  # Based on the epsilon/sigma definition of this potential
  # r float: Interionic separation
  import numpy as np
  lj = (E * (((n/m - n) * (S/r)**m) - ((m/m - n) * (S/r)**n))+(((1*1*e**2)/(4*np.pi*E0*r))*6.242e+18))
  return lj

# Optimisation routine
def brute_force_buck(rho_range, A_range, C_range, r_list)
  # I defined best_area as infinite and best_params as none, so I can use brute-force

  # rho_range tuple: a tuple defining the start, end and spacing (start, end, spacing) of the sample spcase of rho 
  # A_range tuple: a tuple defining the start, end and spacing (start, end, spacing) of the sample spcase of A
  # C_range tuple: a tuple defining the start, end and spacing (start, end, spacing) of the sample spcase of C
  # r_list list: list of the interionic separations to sample over

  from scipy.integrate import simps
  import numpy as np

  best_area = float('inf')
  best_params = None
  for p in np.linspace(rho_range):
      for A in np.linspace(A-range):
          for C in np.linspace(c_range):
           
              # this line of code defines result_list by using A, p and C values and every r values in i_list and uses mm_potential function to fill result_list list
              result_list = buckingham_potential(r, A, rho, C) for r in r_list]


              y_difference = np.abs(np.interp(x_values, i_list, result_list) - (np.interp(x_values, smooth_numbers, smooth_distances)))

              area = simps(y_difference, x=x_values)

              # We are checking if the area we found with the current values of A, p and C is the best one yet. If so, we change the best_area value with area value we found.
              if area < best_area:
                  best_area = area
                  best_params = (A, p, C)

  return best_params
