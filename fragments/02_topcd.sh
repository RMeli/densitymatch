#!/bin/bash

# Compute fragments PCDs

libdir="VEHICLe-good"

n=$(ls ${libdir}/*.sdf | wc -l)

for frag in $(ls VEHICLe-good/*.sdf)
do
    echo ${frag}

    outname=$(echo ${frag} | sed "s/.sdf/.pcd/g")
    time singularity run --nv --app python ../development/densitymatch.sif \
       ../molgrid_to_pcd.py ${frag} --ligmap ../files/ligmap -o ${outname}
done