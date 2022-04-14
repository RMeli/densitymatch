#!/bin/bash

# Generate conformers

#python conformers.py data/DSIpoised.csv
#python conformers.py data/VEHICLe.csv

mkdir -p "VEHICLe-good"
python conformers.py data/VEHICLe_good.csv -d "VEHICLe-good"
