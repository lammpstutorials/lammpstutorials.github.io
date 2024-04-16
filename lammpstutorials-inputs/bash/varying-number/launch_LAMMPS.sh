#!/bin/bash
set -e

lmp="/home/simon/Softwares/lammps-2Aug2023/src/lmp_serial"

for nb2 in 1 9 81 729
do
    echo 'nb2 = '${nb2}
    ${lmp} -in input.lammps -var nb2 ${nb2} -var rdm2 $RANDOM
    folder=nb${nb2}
    mkdir ${folder}
    cp dump.lammpstrj ${folder}
done