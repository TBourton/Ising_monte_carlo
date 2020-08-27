#!/bin/bash --login

#BSUB -J Ising
#BSUB -W 00:10
#BSUB -n 40
#BSUB -x

module purge
module load compiler/gnu-4.8.0
module load mpi/openmpi-1.6.4

x=100
mpic++ -O3 ising_mpi.C random.C -o ising_mpi
for i in `cat nproc.list`; do
        mkdir WS${i}
        let L ="$x*${i}"
	sed "s/DIM/$L/g" ws_input.prov | sed "s/NPROC/${i}/g" > \
		WS${i}/input
	cp ising_mpi WS${i}
	cd WS${i}
        echo Beginning on Nproc = ${i} dim=$L
	mpirun -np ${i} ./ising_mpi < input > out
        echo Nproc = ${i} dim=$L Finished
	cd ..
done		
