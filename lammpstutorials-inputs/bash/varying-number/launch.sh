#!/bin/bash
set -e

lmp="/home/simon/Softwares/lammps-2Aug2023/src/lmp_serial"

for i in 1 9 81 729
do
    ${lmp} -in input.lammps -var nb2 ${i}
    folder=nb${i}
    mkdir ${folder}
    cp dump.lammpstrj ${folder}
done