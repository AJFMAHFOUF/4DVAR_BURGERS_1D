#!/bin/ksh
set -evx
#-------------------------------------------------------------------
# 4D-Var assimilation with 1D Burgers model (advection-diffusion)
#-------------------------------------------------------------------
EXPID=$1
# Execute assimilation
cd ../src
main
# Store results
mv fort.173 ../data/cost_function_gradient_${EXPID}.dat
mv fort.170 ../data/output_24h_${EXPID}.dat
mv fort.171 ../data/output_48h_${EXPID}.dat
