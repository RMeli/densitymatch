"""
Utility functions originally written in Jupyter Notebooks.
"""

# FIXME: Package things properly
# Allow to import scripts in root directory
import sys
sys.path.append("..")

from score_pcd import fit_and_score

import numpy as np

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

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
    coords_aug[:,:3] = coords

    # Compute new (transformed) coordinates
    coords_new = np.matmul(tran, coords_aug.T)[:3,:].T

    # Add new coordinates as conformer
    # Set new conformer index to toConfId
    n_atoms = mol.GetNumAtoms()
    conf = Chem.Conformer(n_atoms)
    conf.SetId(toConfId)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords_new[i,:])

    # Conformer ID manually assigned above, avoid automatic assignement here
    _ = mol.AddConformer(conf, assignId=False)

class AlignShow:
    """
    Align molecules pairs within a list.

    Store scroes and aligned conformers for visualisation (cached).
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
            self._align(idx1, idx2)

        mol1, mol2 = self.mols[idx1], self.mols[idx2]

        # Create view
        p = py3Dmol.view()

        # Select correct conformers
        # mols[idx1] is aligned to mols[idx1]
        p.addModel(Chem.MolToMolBlock(mol1, confId=idx1),'sdf') # mol1 original coordinates
        p.addModel(Chem.MolToMolBlock(mol1, confId=idx2),'sdf') # mol1 aligned to mols[idx2]
        p.addModel(Chem.MolToMolBlock(mol2, confId=idx2),'sdf') # mol2 original coordinates

        p.setStyle({"model": 0}, {'stick':{'colorscheme':'lightgreyCarbon'}})
        p.setStyle({"model": 1}, {'stick':{'colorscheme':'redCarbon'}})
        p.setStyle({"model": 2}, {'stick':{'colorscheme':'greyCarbon'}})

        p.zoomTo()

        return p

def show_molecule_idx(idx, mols):
    """
    Draw a particular molecule from a list.
    """
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mols[idx], confId=0), 'sdf')
    p.setStyle({"model": 0}, {'stick':{'colorscheme':'lightgreyCarbon'}})
    p.zoomTo()
    p.show()


def show_conformer_idx(idx, mol):
    """
    Draw a conformer from a molecule.
    """
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mol, confId=idx), 'sdf')
    p.setStyle({"model": 0}, {'stick':{'colorscheme':'lightgreyCarbon'}})
    p.zoomTo()
    p.show()

def show_all_conformers(mol):
    """
    Draw all conformers for a molecule.
    """
    p = py3Dmol.view()
    for i in range(mol.GetNumConformers()):
        mb = Chem.MolToMolBlock(mol, confId=i)
        p.addModel(mb, 'sdf')
        p.setStyle({"model": i}, {'stick':{'colorscheme':'lightgreyCarbon'}})
    p.zoomTo()
    p.show()