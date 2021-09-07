"""
Fit atoms into atomic densities.

Based on liGAN.

TODO: ADD LICENSE
"""


import molgrid
from rdkit import Chem
from openbabel import pybel

#  FIXME: Make liGAN a proper python package
import sys
sys.path.append("../ligan-EVOTEC") # Evotec version of liGAN (master branch)

# liGAN-related imports
import fitting
from  atom_grids import AtomGrid
from atom_types import get_channels_from_map

if __name__ == "__main__":
    import argparse
    import numpy as np
    import torch

    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    p = argparse.ArgumentParser(
        description="Fit atoms in density difference."
    )
    p.add_argument("numpy", type=str, help="NumPy file")
    p.add_argument(
        "-r", "--resolution", type=float, default=0.5, help="Grid resolution"
    )
    p.add_argument(
        "-m", "--ligmap", type=str, default="files/ligmap", help="Ligand types file"
    )
    p.add_argument(
        "-o", "--output", type=str, default="diff.pcd", help="Output file",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")

    args = p.parse_args()

    # Load grid from NumPy array and change to PyTorch tensor
    npgrid = np.load(args.numpy)
    assert npgrid.shape[0] == 1

    # Load center
    c = np.load(args.numpy.replace(".npy", "-c.npy"))

    if args.verbose:
        print("Grid size:", npgrid.shape)

    #grid = torch.from_numpy(npgrid)
    #grid.to(device)

    # Create GNINA typer
    typer = molgrid.FileMappedGninaTyper(args.ligmap) # Called "lig_map" in liGAN
    lig_channels = get_channels_from_map(typer, name_prefix='Ligand')

    # Move grid to atom grid
    # AtomGrid is a 3D representation of a molecular structure
    # Assumes only one batch with one element (single grid)
    # Need to work with NumPy arrays, not PyTorch tensors
    # Both the values and the center have to be NumPy arrays
    ag = AtomGrid(
            values=npgrid[0], # Use NumPy array to avoid UserWarning
            channels=lig_channels,
            center=c,
            resolution=args.resolution,
        )

    # From the docstring of simple_atom_fit:
    #   Types are ignored
    #   The number of  atoms of each type is inferred from the density
    # Molecule returned is an AtomStruct
    mol, bestgrid = fitting.simple_atom_fit(
        mgrid=ag,
        types=None,
        iters=25,
        tol=0.01,
        grm=1.0,
    )

    print(mol)
    rdmol = mol.to_rd_mol()
    print(rdmol)
    # Get OBMol from AtomStruct
    # Convert OBMol to pybel Molecule for easy output
    obmol = pybel.Molecule(mol.to_ob_mol())
    obmol.write("sdf", "obmol.sdf", overwrite=True)
