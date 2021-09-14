#!/bin/bash

# Generate conformers and compute corresponding PCDs

python conformers.py DSIpoised.csv
python conformers.py VEHICLe.csv

#for frag in $(ls DSIpoised/*.sdf)
#do
#    echo ${frag}
#    outname=$(echo ${frag} | sed "s/.sdf/.pcd/g")
#    time singularity run --nv --app python ../development/densitymatch.sif \
#        ../molgrid_to_pcd.py ${frag} --ligmap ../files/ligmap -o ${outname}
#done