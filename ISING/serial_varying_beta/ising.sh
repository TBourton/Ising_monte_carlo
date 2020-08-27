#!/bin/sh

L=40

 g++ -O3 ising.C random.C -o ising

for i in `cat beta.list`; do
	mkdir B${i}
	sed "s/DIM/$L/g" input.prov | sed "s/BETA/${i}/g" > \
		B${i}/input
	cp ising B${i}
	cd B${i}
	./ising < input > out
	cd ..
done		
