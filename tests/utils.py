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

class AlignShow:

    def __init__(self, mols, pcds):
        """
        Store a series of RDKit molecules and corresponding point clouds.
        """
        assert len(mols) == len(pcds)

        self.mols = mols
        self.pcds = pcds

        # Number of molecules
        self.n = len(mols)

        # Cache for aligned molecules
        self.scores = {}


    def align(self, idx1, idx2):
        """
        Align two molecules based on their PCD representation.
        Uses fit_and_score.
        """
        assert 0 <= idx1 < self.n
        assert 0 <= idx2 < self.n

        if (idx1, idx2) in self.scores:
            # Molecules already aligned
            return self.scores[(idx1, idx2)]

        pcd1, pcd2 = self.pcds[idx1], self.pcds[idx2]
        mol1, mol2 = self.mols[idx1], self.mols[idx2]
    
        fit, cfit, tran = fit_and_score((pcd1, pcd2), voxel_size=0.5, threshold=0.5)

        # Get coordinates to transform
        coords = mol1.GetConformer(0).GetPositions()

        # Augment coordinates with ones
        coords_aug = np.ones((coords.shape[0], 4))
        coords_aug[:,:3] = coords

        # Compute new (transformed) coordinates
        coords_new = np.matmul(tran, coords_aug.T)[:3,:].T

        # Add new coordinates as conformer
        n_atoms = mol1.GetNumAtoms()
        conf = Chem.Conformer(n_atoms)
        for i in range(n_atoms):
            conf.SetAtomPosition(i, coords_new[i,:])

        _ = mol1.AddConformer(conf, assignId=True)

        # Cache scores
        self.scores[(idx1, idx2)] = cfit.fitness

        return self.scores[(idx1, idx2)]

    def show(self, idx1, idx2):
        """
        Show original molecules and aligned structures.
        """
        assert 0 <= idx1 < self.n
        assert 0 <= idx2 < self.n

        if not (idx1, idx2) in self.scores:
            self._align(idx1, idx2)

        mol1, mol2 = self.mols[idx1], self.mols[idx2]

        # Create view
        p = py3Dmol.view()
        p.addModel(Chem.MolToMolBlock(mol1, confId=0),'sdf')
        p.addModel(Chem.MolToMolBlock(mol1, confId=1),'sdf')
        p.addModel(Chem.MolToMolBlock(mol2, confId=0),'sdf')

        p.setStyle({"model": 0}, {'stick':{'colorscheme':'lightgreyCarbon'}})
        p.setStyle({"model": 1}, {'stick':{'colorscheme':'redCarbon'}})
        p.setStyle({"model": 2}, {'stick':{'colorscheme':'greyCarbon'}})

        p.zoomTo()

        return p

def show_molecule_idx(idx, mols):
    p = py3Dmol.view()
    p.addModel(Chem.MolToMolBlock(mols[idx], confId=0), 'sdf')
    p.setStyle({"model": 0}, {'stick':{'colorscheme':'lightgreyCarbon'}})
    p.zoomTo()
    p.show()
