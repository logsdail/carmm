#!/bin/bash

awk 'BEGIN {energy=0; level=0; stack[1]=0; occ=0; stack_occ[1]=0} $1 == "Vector" {energy=$4$5; split(energy,just_energy,"="); occ=$3; split(occ,just_occ,"="); level=$2-100; stack[level]=just_energy[2]; stack_occ[level]=just_occ[2]} END {for (i=1; i<=347; i++) {print stack[i] " " stack_occ[i]}}' nwchem.nwout > wonky.dos

sed -i 's/D/E/g' wonky.dos

awk '{print "(" $1*27.21138602 ", " $2 ")"}' wonky.dos > file.dos

rm wonky.dos
