#!/bin/bash

touch time.txt
touch err.txt


for i in 128 256 512 1024 2048 4096 8192 16384
do 
	sed -n 3p *N$i*.o* >> err.txt
	sed -n 4p *N$i*.o* >> time.txt
done
