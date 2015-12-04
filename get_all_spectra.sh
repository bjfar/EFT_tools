#!/bin/bash

# Runs recoil_generator.m for several isotopes of Xenon

for i in 128 129 130 131 132 134 136; do
    echo "**************************************"
    echo " Generating recoil spectra for Xe"$i"..."
    echo "**************************************"
    ./recoil_generator.m . $i NR
    #./recoil_generator.m . $i use-relativistic
done
