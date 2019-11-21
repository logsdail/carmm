#!/bin/bash

# FILENAME
if [ "${1}" ]; then
input_file=${1}
else
echo "Error1 : No input filename"
exit
fi

# CONVERT ANGSTROM TO BOHR RADIUS
cutoff="0"
bohr_radius="0.52917721092"
if [ "${2}" ]; then
  if [ "${2}" == "norm" ]; then
  normalise="y"
  else
  normalise="n"
    if [[ "${2}" =~ ^[-+]?[0-9]+([.][0-9]+)?$ ]] ; then
      cutoff=`echo ${2}/${bohr_radius} | bc -l`
    fi
  fi
fi

# NORMALISE INPUTS
if [ "${3}" ]; then
  if [ "${3}" == "norm" ]; then
  normalise="y"
  else
  normalise="n"
  fi
fi

if [ -f ${input_file} ]; then 
# LINES IN FILE
lines=`wc -l ${input_file} | awk '{print $1}'`
echo "Input ${input_file} has ${lines} lines"
if [ "${normalise}" == "y" ]; then
echo "Normalising the Gaussian functions"
fi
if ! [ "${cutoff}" == "0" ]; then
echo "Cutoff defined as ${cutoff} Bohr"
fi

# STRINGS OF GNUPLOT INPUT
plot=""
title=""
full_line=""
orbital=""
normaliser=""
# TESTING OF PI IN NORMALISER
# pi=`echo "scale=10; 4*a(1)" | bc -l`
# pi_cubed=`echo "${pi}^3" | bc -l`

# COUNTER
count=1
finish=$lines

while [ ${count} -le ${finish} ]; 
  do
    input_line=`awk NR==${count} ${input_file}`
    input_line_length=`echo ${input_line} | wc -m | awk '{print $1}'`
    first_word=`echo ${input_line} | awk '{print $1}'`
    first_character=${first_word:0:1}
    ((end_marker=`echo ${input_line} | grep -i END | wc -l`))
  
    if [ ${input_line_length} -le 2 ]; then
    ((count=$count+1))
    elif [ "${first_character}" == "#" ]; then    
    ((count=$count+1)) 
    elif [ "${first_word}" == "basis" ]; then
    ((count=$count+1))
#    elif [ ${end_marker} -eq 1 ]; then
#    # INSERT METHOD TO REMOVE FINAL COMMA HERE
#    full_line=`echo ${full_line} | sed '$s/.$//'`
#    ((count=$count+1))
    else
# This just did decimal number. New setup does exponential notation as well
#      if ! [[ "${first_word}" =~ ^[-+]?[0-9]+([.][0-9]+)?$ ]] ; then 
       if ! [[ "${first_word}" =~ ^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
        if [ "${plot}" == "" -a "${title}" == "" ]; then
          full_line=""
        else
          full_line="${full_line}${plot}${title}, "
        fi
        title=" title '${input_line:0:${input_line_length}}'"
        orbital=${first_character}
        plot=""
        ((count=$count+1))
      else
        if [ "${plot}" != "" ]; then
          plot=${plot}" + "
        fi
        prefactor=`echo ${input_line} | awk '{print $1}'`
        exponent=`echo ${input_line} | awk '{print $2}'`

        if [ "${normalise}" == "y" ]; then
          # PUT IN NORMALISING FACTORS HERE
          # TAKEN FROM MODERN QUANTUM CHEMISTRY, SZABO AND OSTLUND
          # SHOULD I MODIFY THESE WITH A SQRT AS SEEN ON GAUSSIAN ORBITAL WIKI?
          if [ "${orbital}" == "S" -o "${orbital}" == "s" ]; then
            normaliser="(((2**0.75)*(${exponent}**0.75))/(pi**0.75))"
          elif [ "${orbital}" == "P" -o "${orbital}" == "p" ]; then
            normaliser="(((2**1.75)*(${exponent}**1.25))/(pi**0.75))"
          elif [ "${orbital}" == "D" -o "${first_character}" == "d" ]; then
            normaliser="(((2**2.75)*(${exponent}**1.75))/(pi**0.75))"
          elif [ "${orbital}" == "F" -o "${first_character}" == "f" ]; then
            normaliser="(((2**3.75)*(${exponent}**2.25))/(pi**0.75))"
          else
            echo "Error2 : Normalisation error. Not designed to deal with ${orbital}-orbitals"
          fi
          equation="(${normaliser}*${prefactor}*exp(-${exponent}*(x**2)))" 
        else
          equation="(${prefactor}*exp(-${exponent}*(x**2)))"
        fi
        plot=${plot}${equation}
        ((count=$count+1))
      fi
    fi
  done

# INSERT METHOD TO REMOVE FINAL COMMA HERE
full_line=`echo ${full_line} | sed '$s/.$//'`
# ORGANISE DISPLAY
plotting_setup="set xlabel 'Bohr Radius'; set ylabel 'Wavefunction';"
if ! [ "${cutoff}" == "0" ]; then
plotting_setup="${plotting_setup} set arrow from ${cutoff},graph(0,0) to ${cutoff},graph(1,1) nohead; \
                set arrow from -${cutoff},graph(0,0) to -${cutoff},graph(1,1) nohead; \
                set xr [((-${cutoff}-5)):((${cutoff}+5))];"
fi
# echo ${plotting_setup}
echo "${plotting_setup} plot ${full_line}" | gnuplot -persist 
# echo ${full_line}
else
echo "Error3 : ${input_file} does not exist"
fi
