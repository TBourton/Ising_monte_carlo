#!/bin/sh

x=40
mpic++ -O3 ising_mpi.C random.C -o ising_mpi
for i in `cat nproc.list`; do
        mkdir WS${i}
        let "L=$x*${i}"
       echo $L
	sed "s/DIM/$L/g" ws_input.prov | sed "s/NPROC/${i}/g" > \
		WS${i}/input
	cp ising_mpi WS${i}
        cp totalmachinefile WS${i}
	cd WS${i}
        echo Beginning on Nproc = ${i} dim=$L
	mpirun -np ${i} -machinefile totalmachinefile ./ising_mpi < input > out
        echo Nproc = ${i} dim=$L Finished
	cd ..
done		
