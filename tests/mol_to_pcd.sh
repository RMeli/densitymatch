#!/bin/bash


for sdf in $(ls ligands/BRD4/*.sdf)
do
    fout="${sdf%.*}.pcd"

    singularity run --nv --app python ../containers/densitymatch.sif \
        ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -r 0.5 -o ${fout} -v
done