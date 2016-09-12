#!/bin/bash

# Runs recoil_generator.m for several isotopes of Xenon

# Location of DMFormFactor package
DMFFdir="/home/farmer/mathematica/DMFormFactor_13086288"

# Output path (where spectrum tables will be written)
outdir="recoil_spectrum_tables"

for isotope in 128 129 130 131 132 134 136; do
#for isotope in 131; do
    #echo "**************************************"
    #echo " Generating recoil spectra for Xe"$isotope"..."
    #echo "**************************************"
    ./recoil_generator.m $DMFFdir $outdir $isotope NR nocheck noHelm 
    echo $?
    if ! echo $?; then 
        echo "Error detected, aborting loop."
        break; 
    fi
    #./recoil_generator.m $DMFFdir $outdir $isotope R nocheck noHelm
    #if ! echo $?; then 
    #    echo "Error detected, aborting loop."
    #    break; 
    #fi
    #./recoil_generator.m $DMFFdir $outdir $isotope NR nocheck UseHelm 
    #echo $?
    #if ! echo $?; then 
    #    echo "Error detected, aborting loop."
    #    break; 
    #fi
    #./recoil_generator.m $DMFFdir $outdir $isotope R nocheck UseHelm
    #if ! echo $?; then 
    #    echo "Error detected, aborting loop."
    #    break; 
    #fi
done
