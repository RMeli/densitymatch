"""
Fit atoms into atomic densities.

Based on liGAN.

TODO: ADD LICENSE
"""


import molgrid
from rdkit import Chem
from rdkit.Chem import rdchem
from openbabel import pybel

#  FIXME: Make liGAN a proper python package
import sys

sys.path.append("../ligan-EVOTEC")  # Evotec version of liGAN (master branch)

# liGAN-related imports
import fitting
from atom_grids import AtomGrid
from atom_types import get_channels_from_map
from molecules import rd_mol_to_ob_mol


def nearest_atoms(rdmol1, rdmol2):
    """
    Given two different molecules, determine the closest atoms.
    """
    mindist2 = np.inf
    midx1, midx2 = -1, -1

    for atom1 in rdmol1.GetAtoms():
        idx1 = atom1.GetIdx()
        pos1 = rdmol1.GetConformer().GetAtomPosition(idx1)

        for atom2 in rdmol2.GetAtoms():
            idx2 = atom2.GetIdx()
            pos2 = rdmol2.GetConformer().GetAtomPosition(idx2)

            d2 = (
                (pos1.x - pos2.x) ** 2 + (pos1.y - pos2.y) ** 2 + (pos1.z - pos2.z) ** 2
            )

            print(np.sqrt(d2))

            if d2 < mindist2:
                mindist2 = d2
                midx1 = idx1
                midx2 = idx2

    return midx1, midx2, np.sqrt(mindist2)


def connectMols(mol1, mol2, idx1, idx2, d):
    """
    Connect to molecules by creating a single bond between given atoms (or one of their
    neighbours if they overlap).

    Notes
    -----
    Original function from:
    https://github.com/molecularsets/moses/blob/master/moses/baselines/combinatorial.py
    """
    atom1 = mol1.GetAtomWithIdx(idx1)
    # atom2 = mol2.GetAtomWithIdx(idx2)

    # Combine the two molecules together and make it editable
    combined = Chem.CombineMols(mol1, mol2)
    emol = Chem.EditableMol(combined)

    if d < 0.25:  # Very short distance, atoms overlap

        # Get nearest neighbour of the overlapping atom in mol1
        # mol1 here is supposed to be the reconstructed side chain
        neighbor1idx = atom1.GetNeighbors()[0].GetIdx()

        # Add bond between nearest neighbour of overlapping atom in mol1
        # And overlapping atom in mol2
        emol.AddBond(
            neighbor1idx, idx2 + mol1.GetNumAtoms(), order=rdchem.BondType.SINGLE
        )

        # Remove overlapping atom in mol1
        emol.RemoveAtom(idx1)
    else:
        # Add bond between nearest atoms in mol1 and mol2
        emol.AddBond(idx1, idx2 + mol1.GetNumAtoms(), order=rdchem.BondType.SINGLE)

    mol = emol.GetMol()

    return mol


if __name__ == "__main__":
    import argparse
    import numpy as np
    import torch

    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    p = argparse.ArgumentParser(description="Fit atoms in density difference.")
    p.add_argument("numpy", type=str, help="NumPy file")
    p.add_argument("scaffold", type=str, help="SDF file")
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

    # grid = torch.from_numpy(npgrid)
    # grid.to(device)

    # Create GNINA typer
    typer = molgrid.FileMappedGninaTyper(args.ligmap)  # Called "lig_map" in liGAN
    lig_channels = get_channels_from_map(typer, name_prefix="Ligand")

    # Move grid to atom grid
    # AtomGrid is a 3D representation of a molecular structure
    # Assumes only one batch with one element (single grid)
    # Need to work with NumPy arrays, not PyTorch tensors
    # Both the values and the center have to be NumPy arrays
    ag = AtomGrid(
        values=npgrid[0],  # Use NumPy array to avoid UserWarning
        channels=lig_channels,
        center=c,
        resolution=args.resolution,
    )

    # From the docstring of simple_atom_fit:
    #   Types are ignored
    #   The number of  atoms of each type is inferred from the density
    # Molecule returned is an AtomStruct
    asmol, bestgrid = fitting.simple_atom_fit(
        mgrid=ag,
        types=None,
        iters=25,
        tol=0.01,
        grm=1.0,  # Gaussian radius multiplier?
    )

    # Get OBMol from AtomStruct
    # Convert OBMol to pybel Molecule for easy output
    # OBMol does not contain bonds at this stage
    obmol = pybel.Molecule(asmol.to_ob_mol())
    obmol.write("sdf", "obmol_fitted.sdf", overwrite=True)

    # Reconstruct OpenBabel molecule from AtomStruct
    obmol, mismatches, visited = fitting.make_obmol(asmol)
    obmol.write("sdf", "obmol_bonds.sdf", overwrite=True)
    print(type(obmol))

    # Convert fitted OpenBabel molecule to RDKit
    # Need to use the OBMol object instead of the pybel molecule
    rdmol = fitting.convert_ob_mol_to_rd_mol(obmol.OBMol)
    rdmol = Chem.RemoveHs(rdmol)  # Remove hydrogens

    # Get original fragment fitted to the density (and removed)
    rdscaffold = next(Chem.SDMolSupplier(args.scaffold, removeHs=True))

    idx1, idx2, d = nearest_atoms(rdmol, rdscaffold)
    print(idx1, idx2, d)

    rdmolfinal = connectMols(rdmol, rdscaffold, idx1, idx2, d)

    # Convert back to OpenBabel, which does not care if things go wrong...
    obmolfinal = pybel.Molecule(rd_mol_to_ob_mol(rdmolfinal))
    obmolfinal.write("sdf", "obmolfinal.sdf", overwrite=True)

    # Try to output RDKkit molecule as well, as MOL file
    # MolToMolBlock allows to use kekulize=False
    with open("rdmolfinal.mol", "w") as fout:
        rdmolfinalh = Chem.AddHs(rdmolfinal, addCoords=True)
        fout.write(Chem.MolToMolBlock(rdmolfinalh, kekulize=False))
