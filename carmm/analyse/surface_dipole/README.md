# surface_dipole
Program to calcule the dipole through a slab. Works with CHGCAR or XYZQ file. Examples provided in "Example" directory.

To compile: 

- type "make"

To use: 

- make sure you have either:
  - a CHGCAR file (VASP output) or 
  - XYZQ file (xyzq.txt, no buffer lines, with N atom lines each containing coordinates x, y, z and charge q. Additionally can includes a uc.txt with unit cell vectors)

- type "./surface_dipole" in working directory

Optional runtime arguments:

--calculate: gives you 2nd and 3rd monents, quadrupole and octopole and eigenvalues of quadrupole (The latter requires lapack to be compiled in - enable LDFLAGS in Makefile)

--compensate: will add a layer of point charges above/below the slab system (z-axis) in order to neutralise a system dipole

--verbose: will give you much more information, including geometric coordinates for compensated systems

--splice: will split the slab in the middle and only perform calculations on the bottom half (if e.g. you have a quadrupole through your slab)

--splice_top: same as above but will just use the top half of the slab

--splice_twice: will splice the system twice, only using the bottom/top quarter of the slab in calculations

Got a question? Email LogsdailA@cardiff.ac.uk
