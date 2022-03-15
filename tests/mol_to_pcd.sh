#!/bin/bash

# ------------------------------------------------------------------------------
# Pre-compute point clouds (PCD) from BRD4 and CDK2 inhibitors using libmolgrid-
# generated densities.
# This is done on translated molecules, so that the point cloud-based alignment
# can be assessed.
# ------------------------------------------------------------------------------

# for ligand in "BRD4" "CDK2"
# do
#     for sdf in $(ls ligands/${ligand}/*.sdf)
#     do
#         fout="${sdf%.*}.pcd"

#         singularity run --nv --app python ../development/densitymatch.sif \
#             ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -r 0.5 -o ${fout} -v
#     done
# done

python transform.py ligands/BRD4/ligand-*.sdf
python transform.py ligands/CDK2/*.sdf

for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ligands/${ligand}/*_tran.sdf)
    do
        fout="${sdf%.*}.pcd"

        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -r 0.5 -o ${fout} -v
    done
done