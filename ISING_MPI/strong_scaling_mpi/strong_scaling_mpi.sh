#!/bin/sh

L=100
mpic++ -O3 ising_mpi.C random.C -o ising_mpi
for i in `cat nproc.list`; do
        mkdir SS${i}
	sed "s/DIM/$L/g" ss_input.prov | sed "s/NPROC/${i}/g" > \
		SS${i}/input
	cp ising_mpi SS${i}
        cp totalmachinefile SS${i}
	cd SS${i}
        echo Beginning on Nproc = ${i}
	mpirun -np ${i} -machinefile totalmachinefile ./ising_mpi < input > out
        echo Nproc = ${i} Finished
	cd ..
done		
