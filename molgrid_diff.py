"""
Compute difference between two molgrids for aligned molecules.
"""

import molgrid
import torch

import numpy as np
import open3d as o3d

from openbabel import pybel

from molgrid_to_pcd import mol_to_grid, grid_to_pcd

if torch.cuda.is_available():
    device = torch.device("cuda")
else:
    device = torch.device("cpu")


def grid_diff(mol1, mol2, dimension, resolution, typer, c=None):
    """
    Compute density difference between two molecules.

    Parameters
    ----------
    mol1:
        Molecule 1
    mol2:
        Molecule 2
    dimension:
        Grid dimension (1D)
    resolution:
        Grid resolution
    typer:
        Atom typer
    c:
        Grid center

    Notes
    -----
    The smaller molecule (less atoms) is always subtracted to the bigger molecule
    (more atoms).
    """
    g1, c1 = mol_to_grid(mol1, dimension, resolution, typer)
    g2, c2 = mol_to_grid(mol2, dimension, resolution, typer, c=c1)

    # FIXME: Is this what we really want?
    # Always subtract the density of the smaller molecule
    # to the density of the bigger molecule
    if len(mol1.atoms) > len(mol2.atoms):
        diff = g1 - g2
    else:
        diff = g2 - g1

    # FIXME: Do we want get the absolute value of the density?
    # Probably not since that would imply re-fitting the removed molecule (known)?
    diff = torch.abs(diff)

    return diff, c1


if __name__ == "__main__":
    import argparse
    import os

    p = argparse.ArgumentParser(
        description="Visualize point cloud (PCD file) with Open3D"
    )
    p.add_argument("sdf", type=str, help="SDF file")
    p.add_argument(
        "-r", "--resolution", type=float, default=0.5, help="Grid resolution"
    )
    p.add_argument(
        "-d", "--dimension", type=float, default=23.5, help="Grid dimension (1D)"
    )
    p.add_argument(
        "-m", "--ligmap", type=str, default="files/ligmap", help="Ligand types file"
    )
    p.add_argument(
        "-o", "--output", type=str, default="diff.pcd", help="Output file",
    )
    p.add_argument("--dx", action="store_true", help="Output grids as DX files")
    p.add_argument(
        "-np", "--numpy", action="store_true", help="Output grid as numpy file"
    )

    # p.add_argument("-v", "--verbose", action="store_true", help="Verbose")

    args = p.parse_args()

    # limbolgrid typer
    typer = molgrid.FileMappedGninaTyper(args.ligmap)

    sdfile = pybel.readfile("sdf", args.sdf)
    obmol1 = next(sdfile)
    obmol2 = next(sdfile)

    diff, c = grid_diff(obmol1, obmol2, args.dimension, args.resolution, typer)

    # Output grids
    if args.dx:
        # https://gnina.github.io/libmolgrid/python/index.html#molgrid.write_dx_grids
        # Grid4f is different from Grid4fCUDA
        # If a function takes Grid4f as input, torch.Tensor need to be moved to the CPU
        print(type(diff), type(typer.get_type_names()), type(c), type(args.resolution))
        molgrid.write_dx_grids(
            prefix="xxx",
            names=typer.get_type_names(),
            grid=diff[0].cpu(),
            center=c,
            resolution=args.resolution,
            scale=1.0,
        )

    if args.numpy:
        np.save(
            args.output.replace(".pcd", ".npy"), diff.cpu().numpy(), allow_pickle=False
        )

    # Convert grid to point cloud
    pcd = grid_to_pcd(diff, c, args.dimension, args.resolution, typer)
    o3d.io.write_point_cloud(args.output, pcd, write_ascii=False)
