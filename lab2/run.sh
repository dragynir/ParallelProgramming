#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=2000m
#PBS -m n
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "Number of MPI process: $MPI_NP"
echo 'File $PBS_NODEFILE:'
cat $PBS_NODEFILE
echo
OMP_NUM_THREADS=1   ./a.out 58 58
OMP_NUM_THREADS=2   ./a.out 58 58
OMP_NUM_THREADS=4   ./a.out 58 58
OMP_NUM_THREADS=8   ./a.out 58 58