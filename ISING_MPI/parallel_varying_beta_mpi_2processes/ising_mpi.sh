#!/bin/sh

L=40

for i in `cat beta.list`; do
	mkdir B${i}
	sed "s/DIM/$L/g" input.prov | sed "s/BETA/${i}/g" > \
		B${i}/input
	cp ising_mpi B${i}
	cd B${i}
	mpirun -np 2 ./ising_mpi < input > out
	cd ..
done		
