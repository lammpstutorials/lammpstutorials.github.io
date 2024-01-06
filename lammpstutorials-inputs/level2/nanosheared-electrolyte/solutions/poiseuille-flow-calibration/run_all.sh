lmp=/work/sgravelle/Softwares/lammps-3Nov2022/src/lmp_mpi

cd force_0.0025/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.005/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.01/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.02/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.04/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

