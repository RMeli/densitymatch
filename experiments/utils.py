"""
Utility functions originally written in Jupyter Notebooks.
"""

# FIXME: Package things properly
# Allow to import scripts in root directory
import sys

sys.path.append("..")

from score_pcd import fit_and_score
from molgrid_to_pcd import sensaas_color_groups_rgb_molgrid

import numpy as np
from numpy.random import default_rng
from scipy.spatial.transform import Rotation as R

# FIXME: Define elsewhere? Allow to control the seed?
rng = default_rng(42)

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

from spyrmsd.rmsd import rmsdwrapper
from spyrmsd.molecule import Molecule


def show(mol1, mol2):
    # Create view
    p = py3Dmol.view()
    # Select correct conformers
    p.addModel(Chem.MolToMolBlock(mol1, confId=0), "sdf")  # mol1 original coordinates
    p.addModel(Chem.MolToMolBlock(mol1, confId=1), "sdf")  # mol1 aligned to mol2
    p.addModel(Chem.MolToMolBlock(mol2, confId=0), "sdf")  # mol2 original coordinates

    p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.setStyle({"model": 1}, {"stick": {"colorscheme": "purpleCarbon"}})
    p.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon"}})

    p.zoomTo()

    return p


def align_and_show(mol1, pcd1, mol2, pcd2):

    rmsd_init = rmsdwrapper(Molecule.from_rdkit(mol1), Molecule.from_rdkit(mol2))

    # Copy original molecule without conformers
    cmol1 = Chem.Mol(mol1, True)

    _, _, tran = fit_and_score((pcd1, pcd2), voxel_size=0.5, threshold=0.5)
    transform_and_add_conformer(mol1, tran, fromConfId=0, toConfId=1)

    cmol1.AddConformer(mol1.GetConformer(1), assignId=True)

    rmsd_final = rmsdwrapper(Molecule.from_rdkit(cmol1), Molecule.from_rdkit(mol2))

    p = show(mol1, mol2)

    return rmsd_init[0], rmsd_final[0], p, mol1


def align(mol1, pcd1, pcd2, hfit=False):
    if hfit:
        gfit, cfit, hhfit, tran = fit_and_score(
            (pcd1, pcd2),
            voxel_size=0.5,
            threshold=0.5,
            color_groups=sensaas_color_groups_rgb_molgrid,
        )
        transform_and_add_conformer(mol1, tran, fromConfId=0, toConfId=1)
        return gfit, cfit, hhfit
    else:
        gfit, cfit, tran = fit_and_score((pcd1, pcd2), voxel_size=0.5, threshold=0.5)
        transform_and_add_conformer(mol1, tran, fromConfId=0, toConfId=1)
        return gfit, cfit

    return gfit, cfit


def show_scaffold(mol, smol):
    # Create view
    p = py3Dmol.view()
    # Select correct conformers
    p.addModel(Chem.MolToMolBlock(mol, confId=0), "sdf")  # mol1 original coordinates
    p.addModel(Chem.MolToMolBlock(smol, confId=1), "sdf")  # mol1 aligned to mol2
    p.addModel(Chem.MolToMolBlock(smol, confId=0), "sdf")  # mol2 original coordinates

    p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.setStyle({"model": 1}, {"stick": {"colorscheme": "purpleCarbon"}})
    p.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon"}})

    p.zoomTo()

    return p


def align_and_show_scaffold(mol, pcd, smol, spcd, sref):
    """
    Align scaffold molecules and show.

    Parameters
    ----------
    mol:
        Full molecule
    pcd:
        Fill molecule point cloud
    smol:
        Murcko scaffold
    spcd:
        Murcko scaffold point cloud
    sref:
        Murcko scaffold reference
    """

    # RMSD between reference and translated Murcko scaffolds
    rmsd_init = rmsdwrapper(Molecule.from_rdkit(smol), Molecule.from_rdkit(sref))

    # Copy original molecule without conformers
    csmol = Chem.Mol(smol, True)

    # Fit scaffold PCD to molecule PCD
    _, _, tran = fit_and_score((spcd, pcd), voxel_size=0.5, threshold=0.5)
    transform_and_add_conformer(smol, tran, fromConfId=0, toConfId=1)

    csmol.AddConformer(smol.GetConformer(1), assignId=True)

    rmsd_final = rmsdwrapper(Molecule.from_rdkit(csmol), Molecule.from_rdkit(sref))

    p = show_scaffold(mol, smol)

    return rmsd_init[0], rmsd_final[0], p, smol


def transform_and_add_conformer(mol, tran, fromConfId=0, toConfId=1) -> None:
    """
    Apply affine transformation to molecule.

    Parameters
    ----------
    mol:
        RDKit molecule
    tran:
        Affine transformation matrix (from Open3D)
    fromConfId:
        Conformer of mol to be transformed (original coordinates)
    toConfId:
        Conformer ID for the transformed conformer (transformed coordinates)
    """
    # Get coordinates to be transformed from conformer
    coords = mol.GetConformer(fromConfId).GetPositions()

    # Augment coordinates with ones
    coords_aug = np.ones((coords.shape[0], 4))
    coords_aug[:, :3] = coords

    # Compute new (transformed) coordinates
    coords_new = np.matmul(tran, coords_aug.T)[:3, :].T

    # Add new coordinates as conformer
    # Set new conformer index to toConfId
    n_atoms = mol.GetNumAtoms()
    conf = Chem.Conformer(n_atoms)
    conf.SetId(toConfId)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords_new[i, :])

    # Ensure that the previous conformer is removed
    # RDKit does not seem to overwrite conformers
    mol.RemoveConformer(toConfId)

    # Conformer ID manually assigned above, avoid automatic assignment here
    _ = mol.AddConformer(conf, assignId=False)


class AlignShow:
    """
    Align molecules pairs within a list.

    Store scores and aligned conformers for visualisation (cached).
    """

    def __init__(self, mols, pcds):
        """
        Store a series of RDKit molecules and corresponding point clouds.

        Parameters
        ----------
        mols:
            List of RDKit molecule
        pcds:
            List of Open3D point clouds (corresponding to mols)
        """
        assert len(mols) == len(pcds)

        self.mols = mols
        self.pcds = pcds

        # Number of molecules
        self.n = len(mols)

        # Cache for aligned molecules
        self.scores = {}

        # Set conformers ID to match molecule id
        for idx, mol in enumerate(self.mols):
            mol.GetConformer().SetId(idx)

        self.all_aligned = False

    def align(self, idx1, idx2):
        """
        Align two molecules based on their PCD representation.
        Uses fit_and_score.

        Parameters
        ----------
        idx1:
            Index of molecule to align
        idx2:
            Index of molecule to which molecule idx1 is being aligned to
        """
        assert 0 <= idx1 < self.n
        assert 0 <= idx2 < self.n

        if (idx1, idx2) in self.scores:
            # Molecules already aligned
            return self.scores[(idx1, idx2)]

        pcd1, pcd2 = self.pcds[idx1], self.pcds[idx2]
        mol1, mol2 = self.mols[idx1], self.mols[idx2]

        fit, cfit, tran = fit_and_score((pcd1, pcd2), voxel_size=0.5, threshold=0.5)

        # Transform coordinates of mol1 (conformer idx1)
        # Conformer index corresponds to the molecule index (modified in constructor)
        # Store new coordinates as conformer with idx2
        transform_and_add_conformer(mol1, tran, fromConfId=idx1, toConfId=idx2)

        # Cache scores
        self.scores[(idx1, idx2)] = cfit.fitness

        return self.scores[(idx1, idx2)]

    def show(self, idx1, idx2):
        """
        Show original molecules and aligned structures.

        Parameters
        ----------
        idx1:
            Index of molecule 1 (show original coordinates and transformed coordinates)
        idx2:
            Index of molecule 2 (show original coordinates)
        """
        assert 0 <= idx1 < self.n
        assert 0 <= idx2 < self.n

        if not (idx1, idx2) in self.scores:
            self.align(idx1, idx2)

        mol1, mol2 = self.mols[idx1], self.mols[idx2]

        # Create view
        p = py3Dmol.view()

        # Select correct conformers
        # mols[idx1] is aligned to mols[idx1]
        p.addModel(
            Chem.MolToMolBlock(mol1, confId=idx1), "sdf"
        )  # mol1 original coordinates
        p.addModel(
            Chem.MolToMolBlock(mol1, confId=idx2), "sdf"
        )  # mol1 aligned to mols[idx2]
        p.addModel(
            Chem.MolToMolBlock(mol2, confId=idx2), "sdf"
        )  # mol2 original coordinates

        p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
        p.setStyle({"model": 1}, {"stick": {"colorscheme": "purpleCarbon"}})
        p.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon"}})

        p.zoomTo()

        return p

    def _align_all(self):
        """
        Align all molecules.
        """
        for i in range(self.n):
            for j in range(self.n):
                self.align(i, j)

        self.all_aligned = True

    def best(self):
        """
        Indices of the best alignment (according to score).
        """
        if not self.all_aligned:
            self._align_all()

        # FIXME: Find better way of eliminating self-alignment
        # The following is faster and more Pythonic, but includes self-alignment
        # max(self.scores, key=self.scores.get)
        best, best_idxs = -1, (-1, -1)
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    continue

                s = self.scores[(i, j)]
                if s > best:
                    best = s
                    best_idxs = (i, j)

        return best, best_idxs

    def best_with(self, idx):
        """
        Best alginment with molecule.
        """
        # TODO" Only align with molecule IDX?
        if not self.all_aligned:
            self._align_all()

        best, best_idxs = -1, (-1, -1)
        for j in range(self.n):
            if j == idx:
                continue

            s = self.scores[(idx, j)]
            if s > best:
                best = s
                best_idxs = (idx, j)

        return best, best_idxs

    def show_best(self):
        """
        Show best alignment (according to score)
        """
        _, (i, j) = self.best()
        return self.show(i, j)

    def save(self, i, j, outfile=None):
        if outfile is None:
            outfile = f"alignment_{i}_{j}.sdf"

        # Make sure conformers to be saved exist
        # Avoid BadConformerId error
        if not (i, j) in self.scores:
            self.align(i, j)

        with Chem.SDWriter(outfile) as w:
            w.write(self.mols[i], confId=j)  # Molecule i aligned to molecule j
            w.write(self.mols[j], confId=j)  # Molecule j


def translate_and_rotate(mol, u=None, t=None, confId=0):
    if u is None:  # Random rotation
        # Random (normalised) axis of rotation
        u = rng.uniform(size=3)
        u /= np.linalg.norm(u)

        # Random angle of rotation
        a = rng.uniform(low=0, high=2 * np.pi)

        # Define rotation
        # Rotation
        r = R.from_rotvec(u * a)
    else:
        # Define rotation from unput u
        r = R.from_rotvec(u)

    if t is None:  # Random translation
        t = rng.uniform(low=-1, high=1, size=3)

    # Affine tranformation
    A = np.zeros((4, 4))
    A[:3, :3] = r.as_matrix()
    A[:3, 3] = t
    A[3, 3] = 1

    # Apply affine transformation to conformer
    rdMolTransforms.TransformConformer(mol.GetConformer(confId), A)


def cut_mol_on_single_bonds(mol):
    """
    Cut molecule over single bonds not in rings and return fragments.

    Original source:
        Project: guacamol_baselines
        Author: BenevolentAI
        File: crossover.py
        License: MIT License
    """
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[*]-;!@[*]")):
        return None

    # List of atom pairs matching SMARTS pattern
    bis = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[*]-;!@[*]")
    )  # single bond not in ring

    fragments = []
    for b in bis:
        # Bond between atom pair
        bs = [mol.GetBondBetweenAtoms(b[0], b[1]).GetIdx()]

        # f = Chem.FragmentOnBonds(mol, bs, addDummies=True, dummyLabels=[(1, 1)])
        f = Chem.FragmentOnBonds(mol, bs, addDummies=False, dummyLabels=[(1, 1)])

        try:
            # Pair of fragments
            mols = Chem.GetMolFrags(f, asMols=True, sanitizeFrags=True)

            # Look for fragments bigger than one heteroatom
            if mols[0].GetNumHeavyAtoms() > 1 and mols[1].GetNumHeavyAtoms() > 1:
                fragments.append(mols)
        except ValueError:
            fragments.appen(None)

    return fragments


def show_molecule_idx(idx, mols):
    """
    Draw a particular molecule from a list.
    """
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mols[idx], confId=0), "sdf")
    p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.zoomTo()
    p.show()


def show_conformer_idx(idx, mol):
    """
    Draw a conformer from a molecule.
    """
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mol, confId=idx), "sdf")
    p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.zoomTo()
    p.show()


def show_all_conformers(mol):
    """
    Draw all conformers for a molecule.
    """
    p = py3Dmol.view()
    for i in range(mol.GetNumConformers()):
        mb = Chem.MolToMolBlock(mol, confId=i)
        p.addModel(mb, "sdf")
        p.setStyle({"model": i}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.zoomTo()
    p.show()


def show_two_mols(mol1, mol2, confId1=0, confId2=0):
    """
    Draw two molecules.
    """
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mol1, confId=confId1), "sdf")
    p.addModel(Chem.MolToMolBlock(mol2, confId=confId2), "sdf")
    p.setStyle({"model": 0}, {"stick": {"colorscheme": "lightgreyCarbon"}})
    p.setStyle({"model": 1}, {"stick": {"colorscheme": "purpleCarbon"}})
    p.zoomTo()
    p.show()
