from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as MS

import argparse
import os

p = argparse.ArgumentParser(description="Optimise and score point clouds")
p.add_argument("mols", type=str, nargs="+", help="Molecules")
p.add_argument("-s", "--suffix", default="murcko", type=str, help="Output suffix")
args = p.parse_args()

for fmol in args.mols:
    # TODO: Use functionailty from utils
    s = Chem.SDMolSupplier(fmol)
    mol = next(s)

    scaffold = MS.GetScaffoldForMol(mol)

    fname, ext = os.path.splitext(fmol)
    outfname = f"{fname}_{args.suffix}{ext}"
    with Chem.SDWriter(outfname) as w:
        w.write(scaffold, confId=0)
