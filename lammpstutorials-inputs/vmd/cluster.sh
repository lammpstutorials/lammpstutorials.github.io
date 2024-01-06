#!/bin/bash
#OAR -n nemd
#OAR -l /nodes=1/cpu=1/core=4,walltime=10:00:00
#OAR --stdout log.out
#OAR --stderr log.err
#OAR --project tamtam

lmp=/home/gravells/softwares/lammps-23Jun2022/src/lmp_mpi

mpirun -np 4 ${lmp} -in equilibrate.lammps
mpirun -np 4 ${lmp} -in mix.lammps
mpirun -np 4 ${lmp} -in video.lammps
