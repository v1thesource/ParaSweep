#!/bin/bash

# This script runs the Brouwer diagram script at param value $param and the smallest possible $Convergence.
# It will append the Brouwer diagram data to datadump.txt, and then increment $param by $param_increment

param_name='Temperature'
param_value='500'
param_maximum='2500'
Convergence=1e-40
param_increment='5'

while [ $( bc <<< "$param_value <= $param_maximum" ) == "1" ]
do
  perl defect_analysis3.pl intrinsic > failure.test
  echo $param_name = $param_value
  echo Convergence = $Convergence
  if grep -Fxq "Could not deteremine the Fermi level that gives charge neutrality, recommned you examine your DFT energies" ./failure.test
    then
    sed -i '' "s/Convergence : $Convergence/Convergence : $(python<<<"print($Convergence*10)")/g" intrinsic.input
    Convergence=$(python<<<"print($Convergence*10)")
  else
    cat intrinsic.res >> datadump.txt
    sed -i '' "s/$param_name : $param_value/$param_name : $( bc <<< "$param_value + $param_increment" )/g" intrinsic.input
    param_value=$( bc <<< "$param_value + $param_increment" )
  fi
done
