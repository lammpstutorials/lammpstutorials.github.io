lmp=/work/sgravelle/Softwares/lammps-3Nov2022/src/lmp_mpi

cd force_0.01/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.02/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.03/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.04/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.05/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.06/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.07/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..

cd force_0.08/
    mpirun -np 8 ${lmp} -in input.lammps
cd ..
