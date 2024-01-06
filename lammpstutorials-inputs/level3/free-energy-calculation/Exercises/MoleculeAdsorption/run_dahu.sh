#!/bin/bash
#OAR -n emd
#OAR -l /nodes=1/cpu=1/core=2,walltime=48:00:00
#OAR --stdout log-water.out
#OAR --stderr log-water.err
#OAR --project tamtam

lmp=/home/gravells/softwares/lammps-2Aug2023/src/lmp_mpi

mpirun -np 2 ${lmp} -in input.lammps
