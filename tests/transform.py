from rdkit.Chem import AllChem as Chem

import numpy as np
from scipy.spatial.transform import Rotation as R

import argparse
import os

from utils import transform_and_add_conformer

p = argparse.ArgumentParser(description="Optimise and score point clouds")
p.add_argument("mols", type=str, nargs="+", help="Molecules")
p.add_argument("-s", "--suffix", default="tran", type=str, help="Output suffix")
args = p.parse_args()


rng = np.random.default_rng(seed=42)

for fmol in args.mols:
    s = Chem.SDMolSupplier(fmol)
    mol = next(s)

    # Translation vector
    t = rng.random((3,)) * 10

    # Rotation matrix
    angle = np.pi * rng.random()
    e = rng.random((3,))
    r = R.from_rotvec(angle * e)

    # Define affine transformation matrix
    A = np.zeros((4, 4))
    A[:3, :3] = r.as_matrix()
    A[:3, 3] = t
    A[3, 3] = 1

    transform_and_add_conformer(mol, A)

    fname, ext = os.path.splitext(fmol)
    outfname = f"{fname}_{args.suffix}{ext}"
    with Chem.SDWriter(outfname) as w:
        w.write(mol, confId=1)
