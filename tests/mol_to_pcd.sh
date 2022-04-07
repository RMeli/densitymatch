#!/bin/bash

# ------------------------------------------------------------------------------
# Pre-compute point clouds (PCD) from BRD4 and CDK2 inhibitors using libmolgrid-
# generated densities.
# This is done on translated molecules, so that the point cloud-based alignment
# can be assessed.
# ------------------------------------------------------------------------------

for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ligands/${ligand}/*.sdf)
    do
        fout="${sdf%.*}.pcd"

        # libmolgrid color scheme
        fout="${sdf%.*}_molgrid.pcd"
        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -o ${fout} -v

        # SENSAAS color scheme
        fout="${sdf%.*}_sensaas.pcd"
        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -o ${fout} -v --sensaas
    done
done

python transform.py ligands/BRD4/ligand-*.sdf
python transform.py ligands/CDK2/*.sdf

for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ligands/${ligand}/*_tran.sdf)
    do
        # libmolgrid color scheme
        fout="${sdf%.*}_molgrid.pcd"
        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -o ${fout} -v

        # SENSAAS color scheme
        fout="${sdf%.*}_sensaas.pcd"
        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -o ${fout} -v --sensaas
    done
done
