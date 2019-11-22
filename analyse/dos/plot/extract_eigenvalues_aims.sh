#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Two inputs needed: input file and output file"
  echo "Example usage: ./script input1 output1"
fi

if [ ! -f $1 ]; then
  echo "File: $1 does not exist"
  exit
fi

#if [ -f $2 ]; then
#  echo "File: $2 does exist - remove it!"
#  exit
#fi

input_file=$1
output_file=$2
number_of_states=`grep "| Number of Kohn-Sham states (occupied + empty):" $input_file | awk '{ print $9 }'`
spin_polarised=`grep "Spin-up eigenvalues:" $input_file | wc -l`
if [ $(($spin_polarised)) -gt 0 ]; then
  grep "Spin-up eigenvalues:" $input_file -A $(($number_of_states+3)) | tail -$number_of_states | awk '{ print $4" "$2 }' > $output_file
  grep "Spin-down eigenvalues:" $input_file -A $(($number_of_states+3)) | tail -$number_of_states | awk '{ print $4" "$2 }' >> $output_file
else
  grep "Writing Kohn-Sham eigenvalues." $input_file -A $(($number_of_states+2)) | tail -$number_of_states | awk '{ print $4" "$2 }' > $output_file
fi

echo "Data written to $output_file"

