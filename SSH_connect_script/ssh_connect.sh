#!/bin/sh

for i in `cat nodelist.list`; do
       echo "yes\n" | ssh ${i}
done		
