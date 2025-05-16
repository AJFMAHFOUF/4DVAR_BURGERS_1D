#!/bin/ksh
set -evx
#-------------------------------------------------------------------
# 3D-Var assimilation with 1D Burgers model (advection-diffusion)
#-------------------------------------------------------------------
EXPID=$1
# Execute assimilation
cd ../src
main
# Store results
split -21 fort.173
mv xaa ../data/cost_function_gradient_${EXPID}_00.dat
mv xab ../data/cost_function_gradient_${EXPID}_06.dat
mv xac ../data/cost_function_gradient_${EXPID}_12.dat
mv xad ../data/cost_function_gradient_${EXPID}_18.dat
mv xae ../data/cost_function_gradient_${EXPID}_24.dat
mv fort.170 ../data/output_24h_${EXPID}.dat
mv fort.171 ../data/output_48h_${EXPID}.dat
mv fort.172 ../data/output_00h_${EXPID}.dat
