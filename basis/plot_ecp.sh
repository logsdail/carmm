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
scalar="n"
bohr_radius="0.52917721092"
if [ "${2}" ]; then
  if [ "${2}" == "scale" ]; then
    scalar="y"
  elif [[ "${2}" =~ ^[-+]?[0-9]+([.][0-9]+)?$ ]] ; then
      cutoff=`echo ${2}/${bohr_radius} | bc -l`
  fi
fi

if [ "${3}" -a "${3}" == "scale" ]; then
  scalar="y"
fi

if [ -f ${input_file} ]; then 
# LINES IN FILE
lines=`wc -l ${input_file} | awk '{print $1}'`
echo "Input ${input_file} has ${lines} lines"
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
((contracted_counter=0))

while [ ${count} -le ${finish} ]; 
  do
    input_line=`awk NR==${count} ${input_file}`
    input_line_length=`echo ${input_line} | wc -m | awk '{print $1}'`
    first_word=`echo ${input_line} | awk '{print $1}'`
    first_character=${first_word:0:1}

#    echo ${input_line} ${input_line_length} ${first_word} ${first_character}
  
    if [ ${input_line_length} -le 2 ]; then
    ((count=$count+1))
      if [ "${plot}" == "" -a "${title}" == "" ]; then
        full_line=""
# WE HAVE ALREADY RESET PLOT IF IT IS BLANK
      elif [ "${plot}" != "" ]; then
#        title=" title '${pretitle} ${plot}'"
          title=" title '${pretitle} ${orbital}'"
          full_line="${full_line}${plot}${title}, "
          plot=""
        
          # Increase orbital character
          if [ "${orbital}" == "ul" ]; then
            orbital="s"
          elif [ "${orbital}" == "s" ]; then
            orbital="p"
          elif [ "${orbital}" == "p" ]; then
            orbital="d"
          elif [ "${orbital}" == "d" ]; then
            orbital="f"
          else
            echo "Error2 : Orbital error. Not designed to deal with ECPs beyond ${orbital}-orbitals"
          fi
        fi
    elif [ "${first_character}" == "#" ]; then   
    ((count=$count+1))
    elif [ "${first_word}" == "pseudo" ]; then
    # THIS SHOULD ONLY BE AT THE TOP OF THE FILE
    ((count=$count+1)) 
    elif [ `echo ${first_word} | tr [:upper:] [:lower:]` == "cards" ]; then
    # HACK IF THERE IS NO SPACES
    if [ "${plot}" != "" ]; then
      title=" title '${pretitle} ${orbital}'"
      full_line="${full_line}${plot}${title}, "
      plot=""
    fi
    pretitle=`echo ${input_line} | awk '{print $2}'`
    orbital="ul"
    ((count=$count+1))
    else
# Read in Gaussian ECP
      rexponent=`echo ${input_line} | awk '{print $1}'`
      prefactor=`echo ${input_line} | awk '{print $2}'`
      exponent=`echo ${input_line} | awk '{print $3}'`
      # Check we have a valid Gaussin and not filler
      if [ "${exponent}" != "" ]; then
        if [ "${plot}" != "" ]; then
          plot=${plot}" + "
        fi
        if [ "${scalar}" == "y" ]; then
          # http://www.nwchem-sw.org/index.php/ECP
          equation="((${prefactor}*(x**(${rexponent}-2))*exp(-${exponent}*(x**2)))/(x**2))"
        else
          equation="(${prefactor}*(x**(${rexponent}-2))*exp(-${exponent}*(x**2)))"
        fi
        plot=${plot}${equation}
      fi
      ((count=$count+1))
      # echo ${pretitle} ${input_line}
    fi
  done

# HACK IF THERE IS NO SPACE AT END OF FILE
if [ "${plot}" != "" ]; then
  title=" title '${pretitle} ${count}'"
  full_line="${full_line}${plot}${title}, "
  plot=""
fi

# INSERT METHOD TO REMOVE FINAL COMMA HERE
full_line=`echo ${full_line} | sed '$s/.$//'`
# ORGANISE DISPLAY
plotting_setup="set xlabel 'Bohr Radius'; set ylabel 'Vecp';"
if ! [ "${cutoff}" == "0" ]; then
plotting_setup="${plotting_setup} set arrow from ${cutoff},graph(0,0) to ${cutoff},graph(1,1) nohead; \
                set arrow from -${cutoff},graph(0,0) to -${cutoff},graph(1,1) nohead; \
                set xr [((-${cutoff}-5)):((${cutoff}+5))];"
fi
# echo ${plotting_setup}
echo "${plotting_setup} plot ${full_line}" | gnuplot -persist 
echo ${full_line}
else
echo "Error3 : ${input_file} does not exist"
fi
