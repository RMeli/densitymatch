"""
Fit atoms into atomic densities.

Based on liGAN.

TODO: ADD LICENSE
"""

import numpy as np


import molgrid

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize

from openbabel import pybel

#  FIXME: Make liGAN a proper python package
import sys

sys.path.append("../ligan-EVOTEC")  # Evotec version of liGAN (master branch)

# liGAN-related imports
import fitting
from atom_grids import AtomGrid
from atom_types import get_channels_from_map


def atom_distances(rdmol1, rdmol2):
    """
    Compute distances between all atom pairs of two different molecules.
    """
    distances = []

    # TODO: Use getPositions() to speed this up
    for atom1 in rdmol1.GetAtoms():
        idx1 = atom1.GetIdx()
        pos1 = rdmol1.GetConformer().GetAtomPosition(idx1)

        for atom2 in rdmol2.GetAtoms():
            idx2 = atom2.GetIdx()
            pos2 = rdmol2.GetConformer().GetAtomPosition(idx2)

            d2 = (
                (pos1.x - pos2.x) ** 2 + (pos1.y - pos2.y) ** 2 + (pos1.z - pos2.z) ** 2
            )

            distances.append((idx1, idx2, np.sqrt(d2)))

    distances.sort(key=lambda t: t[-1])

    return distances


def nearest_atoms(rdmol1, rdmol2):
    """
    Given two different molecules, determine the closest atoms.
    """
    adists = atom_distances(rdmol1, rdmol2)

    return adists[0]


def connectMols(mol1, mol2, idx1, idx2, d):
    """
    Connect to molecules by creating a single bond between given atoms (or one of their
    neighbours if they overlap).

    Notes
    -----
    https://colab.research.google.com/github/pschmidtke/blog/blob/master/_notebooks/2021-01-23-grafting-fragments.ipynb#scrollTo=vhZlPkDHILQ7

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


def remove_overlapping(mol1, mol2, threshold=0.3):
    """
    Remove atoms from mol2 overlapping with mol1.
    """

    emol = Chem.EditableMol(mol2)

    # TODO: Optimize
    adists = atom_distances(mol1, mol2)

    to_remove = []
    for _, j, d in adists:
        if d < threshold:
            to_remove.append(j)

    # When atom is deleted the index is lost
    # https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/2017062518163195443411%40gmail.com/
    for r in sorted(to_remove, reverse=True):
        emol.RemoveAtom(r)

    return emol.GetMol()


def fit_atoms(diff, center, resolution, ligmap, verbose=False):
    """
    Fit libmolgrid density obtained as difference between a scaffold/fragment and
    a (generated) density from a larger molecule.
    """
    # Create GNINA typer
    typer = molgrid.FileMappedGninaTyper(ligmap)  # Called "lig_map" in liGAN
    lig_channels = get_channels_from_map(typer, name_prefix="Ligand")

    # If center is a molgrid vector, convert to numpy array
    if isinstance(center, molgrid.float3):
        center = np.array([*center])

    # Move grid to atom grid
    # AtomGrid is a 3D representation of a molecular structure
    # Assumes only one batch with one element (single grid)
    # Need to work with NumPy arrays, not PyTorch tensors
    # Both the values and the center have to be NumPy arrays
    ag = AtomGrid(
        values=diff[0],  # Use NumPy array to avoid UserWarning
        channels=lig_channels,
        center=center,
        resolution=resolution,
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
    if verbose:
        obmol.write("sdf", "atoms.sdf", overwrite=True)

    # Reconstruct OpenBabel molecule from AtomStruct
    obmol, mismatches, visited = fitting.make_obmol(asmol)
    if verbose:
        obmol.write("sdf", "obmol.sdf", overwrite=True)

    # Convert fitted OpenBabel molecule to RDKit
    # Need to use the OBMol object instead of the pybel molecule
    rdmol = fitting.convert_ob_mol_to_rd_mol(obmol.OBMol)
    rdmol = Chem.RemoveHs(rdmol)  # Remove hydrogens

    if verbose:
        # Output final reconstructed molecule
        with open("rdmol.mol", "w") as fout:
            fout.write(Chem.MolToMolBlock(rdmol))

    return rdmol


def molgrid_diff_to_mol(diff, center, resolution, ligmap, scaffold=None, verbose=False):
    """
    Fit atoms to density difference and (eventually) link such fragment to a given
    scaffold.

    Parameters
    ----------
    diff: np.ndarray
        Atomic density (obtained as difference)
    center:
        Grid centers
    resolution:
        Grid resolutioon
    ligmap:
        Path to ligand mapping file
    scaffold:
        Optional scaffold to link to the fitted atoms (reconstructed fragment)
    verbose:
        Flag to print out intermediate states
    """
    # Molecule (fragment) built by fitting atoms into the density difference
    rdmol = fit_atoms(diff, center, resolution, ligmap, verbose)

    if scaffold is not None:
        rdmolfinal = scaffold  # Initialise molecule with scaffold

        # Remove fitted atoms overlapping with scaffold
        # This should create disconnected component
        # Avoid problems with two ring substituents bound together
        rdmol = remove_overlapping(rdmolfinal, rdmol)
        if verbose:
            # Output final reconstructed molecule
            with open("rdmol-noverlap.mol", "w") as fout:
                fout.write(Chem.MolToMolBlock(rdmol))

        # Fitted atoms can be in different part of the scaffold
        # They constitute different disconnected components (fragments)
        for frag in Chem.GetMolFrags(rdmol, asMols=True, sanitizeFrags=True):
            # Get all distance between fragment and scaffold (sorted)
            idx1, idx2, d = nearest_atoms(frag, rdmolfinal)

            # Connect fragment to elaborated scaffold
            rdmolfinal = connectMols(frag, rdmolfinal, idx1, idx2, d)

    else:
        rdmolfinal = rdmol

    # Try to sanitize and normalise molecule
    # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
    rdmolfinal.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(
        rdmolfinal,
        sanitizeOps=(
            Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES
        ),
    )
    rdmolfinal = rdMolStandardize.Normalize(rdmolfinal)

    return rdmolfinal


if __name__ == "__main__":
    import argparse
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
        "-o", "--output", type=str, default="rdmolfinal.mol", help="Output file",
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
        print("Center:", c)

    # Get original fragment fitted to the density (and removed)
    rdscaffold = next(Chem.SDMolSupplier(args.scaffold, removeHs=True))

    # Fit atoms into density difference
    # Link nearest atom from the fit to the scaffold to build whole molecule
    rdmolfinal = molgrid_diff_to_mol(
        npgrid, c, args.resolution, args.ligmap, rdscaffold, args.verbose
    )

    # Output final reconstructed molecule
    with open(args.output, "w") as fout:
        # rdmolfinalh = Chem.AddHs(rdmolfinal, addCoords=True)
        fout.write(Chem.MolToMolBlock(rdmolfinal))
