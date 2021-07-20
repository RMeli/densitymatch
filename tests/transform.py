from rdkit.Chem import AllChem as Chem

import numpy as np
from scipy.spatial.transform import Rotation as R

import argparse
import os

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
    A = np.zeros((4,4))
    A[:3,:3] = r.as_matrix()
    A[:3,3] = t
    A[3,3] = 1
    #print(A)

    # Get coordinates to transform
    coords = mol.GetConformer(0).GetPositions()

    # Augment coordinates with ones
    coords_aug = np.ones((coords.shape[0], 4))
    coords_aug[:,:3] = coords

    # Compute new (transformed) coordinates
    coords_new = np.matmul(A, coords_aug.T)[:3,:].T

    # Add new coordinates as conformer
    n_atoms = mol.GetNumAtoms()
    conf = Chem.Conformer(n_atoms)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords_new[i,:])

    # Add conformer
    _ = mol.AddConformer(conf, assignId=True)

    fname, ext = os.path.splitext(fmol)
    outfname = f"{fname}_{args.suffix}{ext}"
    with Chem.SDWriter(outfname) as w:
        w.write(mol, confId=1)