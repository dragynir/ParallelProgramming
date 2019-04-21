#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=8:mpiprocs=8:mem=2000m
#PBS -m n
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "Number of MPI process: $MPI_NP"
echo 'File $PBS_NODEFILE:'
cat $PBS_NODEFILE
echo
mpirun -machinefile $PBS_NODEFILE -np 1  ./a.out 512 384 256
mpirun -machinefile $PBS_NODEFILE -np 2  ./a.out 512 384 256
mpirun -machinefile $PBS_NODEFILE -np 4  ./a.out 512 384 256
mpirun -machinefile $PBS_NODEFILE -np 8  ./a.out 512 384 256
mpirun -machinefile $PBS_NODEFILE -np 16 ./a.out 512 384 256