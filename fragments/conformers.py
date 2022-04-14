"""
Generate one conformer per fragment, for all fragments in a given library (CSV file)
"""

import pandas as pd

import tqdm

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdDistGeom

if __name__ == "__main__":
    import argparse
    import os

    p = argparse.ArgumentParser(
        description="Create conformers for fragments (SDF and PCD)."
    )
    p.add_argument("library", type=str, help="Fragment library.")
    p.add_argument(
        "-c", "--conformers", type=int, default=1, help="Number of conformers."
    )
    p.add_argument(
        "-d", "--directory", type=str, default=None, help="Output directory."
    )
    p.add_argument("-s", "--seed", type=int, default=42, help="RNG seed.")
    args = p.parse_args()

    libname, _ = os.path.splitext(args.library)

    df = pd.read_csv(args.library)

    # Determine if SMILES columns is called "smiles" or "SMILES"
    for strf in [str.upper, str.lower]:
        smi_col = strf("smiles")
        if smi_col in df.columns:
            break

    PandasTools.AddMoleculeColumnToFrame(df, smi_col, "mol")

    outpath = libname if args.directory is None else args.directory
    print(outpath)

    for idx, row in tqdm.tqdm(df.iterrows(), desc=libname, total=len(df)):
        molh = Chem.AddHs(row.mol, addCoords=True)

        param = rdDistGeom.ETKDGv3()
        param.random_seed = args.seed
        cids = rdDistGeom.EmbedMultipleConfs(molh, args.conformers, param)

        outsdf = os.path.join(outpath, f"fragment_{idx}.sdf")
        w = Chem.SDWriter(outsdf)
        for cid in cids:
            molh.SetProp("CID", str(cid))  # Add conformer ID as molecular property
            w.write(molh, confId=cid)
        w.close()
