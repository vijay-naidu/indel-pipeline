#!/bin/bash

names='ERR502366'
first='_1'
second='_2'
for name in $names
do
echo $name
trim_galore --length 70 --paired -o trim $name$first'.fastq' $name$second'.fastq'
done
echo All done