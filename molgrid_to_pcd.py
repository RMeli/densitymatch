"""
Convert molecule into point cloud using an isosurface of libmolgrid-computed
atomic densities.
"""

import molgrid
import torch

import numpy as np
import open3d as o3d

from openbabel import pybel
from collections import defaultdict

if torch.cuda.is_available():
    device = torch.device("cuda")
else:
    device = torch.device("cpu")

name_to_rgb_sensaas = defaultdict(
    lambda: (0, 0, 255),  # Default value (everything else)
    {
        "AliphaticCarbonXSHydrophobe": (0, 255, 0),
        "AliphaticCarbonXSNonHydrophobe": (0, 255, 0),
        "AromaticCarbonXSHydrophobe": (0, 255, 0),
        "AromaticCarbonXSNonHydrophobe": (0, 255, 0),
        "Bromine": (230, 230, 230),
        "Chlorine": (230, 230, 230),
        "Fluorine": (255, 0, 0),
        "Nitrogen": (255, 0, 0),
        "NitrogenXSAcceptor": (255, 0, 0),
        "NitrogenXSDonor": (255, 0, 0),
        "NitrogenXSDonorAcceptor": (255, 0, 0),
        "Oxygen": (255, 0, 0),
        "OxygenXSAcceptor": (255, 0, 0),
        "OxygenXSDonorAcceptor": (255, 0, 0),
        "Phosphorus": (0, 255, 0),
        "Sulfur": (255, 0, 0),
        "SulfurAcceptor": (255, 0, 0),
        "Iodine": (230, 230, 230),
        "Boron": (0, 255, 0),
    },
)

# Some elements are combines because their propery can change 
# in fragments (change of substituents).
name_to_rgb_molgrid = {
    "AliphaticCarbonXSHydrophobe": (96, 96, 96),
    "AliphaticCarbonXSNonHydrophobe": (96, 96, 96),
    "AromaticCarbonXSHydrophobe": (192, 192, 192),
    "AromaticCarbonXSNonHydrophobe": (192, 192, 192),
    "Bromine": (102, 0, 51),
    "Chlorine": (0, 204, 0),
    "Fluorine": (51, 255, 255),
    "Nitrogen": (0, 0, 255),
    "NitrogenXSAcceptor": (0, 0, 255),
    "NitrogenXSDonor": (0, 0, 255),
    "NitrogenXSDonorAcceptor": (0, 0, 255),
    "Oxygen": (153, 0, 0),
    "OxygenXSAcceptor": (153, 0, 0),
    "OxygenXSDonorAcceptor": (153, 0, 0),
    "Phosphorus": (153, 76, 0),
    "Sulfur": (153, 153, 0),
    "SulfurAcceptor": (153, 153, 0),
    "Iodine": (0, 204, 0),
    "Boron": (153, 76, 0),
}


def mol_to_grid(obmol, dimension, resolution, typer, c=None):
    gm = molgrid.GridMaker(resolution=resolution, dimension=dimension)

    # Grid dimensions (including types)
    gdims = gm.grid_dimensions(typer.num_types())

    # Pre-allocate grid
    # Only one example (batch size is 1)
    grid = torch.zeros(1, *gdims, dtype=torch.float32, device=device)

    # Add hydrogens to molecule
    obmol.addh()

    # This can cause a MemoryError
    # See https://github.com/gnina/libmolgrid/issues/62
    cs = molgrid.CoordinateSet(obmol, typer)

    # Create example
    ex = molgrid.Example()
    ex.coord_sets.append(cs)

    # Compute center
    if c is None:
        c = ex.coord_sets[0].center()  # Only one coordinate set

    # https://gnina.github.io/libmolgrid/python/index.html#the-transform-class
    # Define transformation to fix center
    transform = molgrid.Transform(c, random_translate=0.0, random_rotation=False)

    # Compute grid
    gm.forward(ex, transform, grid[0])

    return grid, c


def _grid_lims(o, L):
    """
    Compute grid limits given grid length and grid origin.
    """
    return o - L / 2.0, o + L / 2.0


def grid_to_pcd(
    grid,
    center,
    dimension,
    resolution,
    typer,
    color_map=name_to_rgb_molgrid,
    surface=True,
):

    if surface:
        # Compute isosurface
        cloud = torch.logical_and(grid[0] >= 0.4, grid[0] <= 0.6)
    else:
        # Use dense representation
        cloud = grid[0] > 0.5

    # Number of voxels in one direction
    n_steps = round(dimension / resolution) + 1

    x1, x2 = _grid_lims(center[0], dimension)
    x = torch.linspace(x1, x2, steps=n_steps)

    y1, y2 = _grid_lims(center[1], dimension)
    y = torch.linspace(y1, y2, steps=n_steps)

    z1, z2 = _grid_lims(center[2], dimension)
    z = torch.linspace(z1, z2, steps=n_steps)

    X, Y, Z = torch.meshgrid(x, y, z)

    # Total number of points in the cloud
    n_points = torch.sum(cloud)

    XYZ = torch.empty(n_points, 3)
    RGB = torch.empty(n_points, 3)

    start = 0
    for i, name in enumerate(typer.get_type_names()):
        idx = cloud[i]  # Indices of cloud points
        rgb = color_map[name]

        # print(f"{name}: {round(torch.sum(idx).item())}")

        step = torch.sum(idx).item()
        stop = start + step

        # print(start, stop, step)
        XYZ[start:stop, 0] = X[idx].reshape((-1,))
        XYZ[start:stop, 1] = Y[idx].reshape((-1,))
        XYZ[start:stop, 2] = Z[idx].reshape((-1,))
        RGB[start:stop, 0] = rgb[0]
        RGB[start:stop, 1] = rgb[1]
        RGB[start:stop, 2] = rgb[2]

        start = stop

    # Scale RGB values
    # RGB values needs to be in [0,1]
    RGB = RGB / 255

    # Manually build point cloud object
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(XYZ.cpu().numpy())
    pcd.colors = o3d.utility.Vector3dVector(RGB.cpu().numpy())

    return pcd


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
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output file (.pcd, .xyzrgb, ...)",
    )
    p.add_argument("--dx", action="store_true", help="Output grids as DX files")
    p.add_argument("--ascii", action="store_true", help="Output file in ASCII format")
    p.add_argument("--sensaas", action="store_true", help="Use SENSAAS color scheme")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")

    args = p.parse_args()

    system = os.path.splitext(os.path.basename(args.sdf))[0]

    if args.output is None:
        args.output = f"{system}.pcd"

    obmol = next(pybel.readfile("sdf", args.sdf))

    if args.verbose:
        print("molecule:\n\t", obmol, end="")

    # limbolgrid typer
    typer = molgrid.FileMappedGninaTyper(args.ligmap)

    if args.verbose:
        print("types:")
        for t in typer.get_type_names():
            print(f"\t{t}")

    # Convert molecule into grid
    grid, center = mol_to_grid(obmol, args.dimension, args.resolution, typer)

    if args.verbose:
        print("grid center:\n\t", tuple(center))
        print("grid shape:\n\t", grid.shape)

    # Output grids
    if args.dx:
        # https://gnina.github.io/libmolgrid/python/index.html#molgrid.write_dx_grids
        # Grid4f is different from Grid4fCUDA
        # If a function takes Grid4f as input, torch.Tensor need to be moved to the CPU
        molgrid.write_dx_grids(
            prefix=f"grids/{system}",
            names=typer.get_type_names(),
            grid=grid[0].cpu(),
            center=center,
            resolution=args.resolution,
            scale=1.0,
        )

    # Convert grid to point cloud
    pcd = grid_to_pcd(
        grid,
        center,
        args.dimension,
        args.resolution,
        typer,
        color_map=name_to_rgb_sensaas if args.sensaas else name_to_rgb_molgrid,
    )

    o3d.io.write_point_cloud(args.output, pcd, write_ascii=args.ascii)
