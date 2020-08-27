#!/bin/bash --login

#BSUB -J Ising
#BSUB -W 01:00
#BSUB -n 40
#BSUB -x

module purge
module load compiler/gnu-4.8.0
module load mpi/openmpi-1.6.4

L=480
mpic++ -O3 ising_mpi.C random.C -o ising_mpi
for i in `cat nproc.list`; do
        mkdir SS${i}
	sed "s/DIM/$L/g" ss_input.prov | sed "s/NPROC/${i}/g" > \
		SS${i}/input
	cp ising_mpi SS${i}
	cd SS${i}
        echo Beginning on Nproc = ${i}
	mpirun -np ${i} ./ising_mpi < input > out
        echo Nproc = ${i} Finished
	cd ..
done		
