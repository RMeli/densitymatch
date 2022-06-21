import sys

root = "../../"
sys.path.append(root)

import os
from openbabel import pybel
from molgrid_to_pcd import mol_to_grid
import molgrid
import argparse

p = argparse.ArgumentParser(
        description="Visualize point cloud (PCD file) with Open3D"
    )
p.add_argument("mol_sdf", type=str, help="Molecule SDF file")
p.add_argument("scaffold_sdf", type=str, help="Scaffold SDF file")
p.add_argument(
    "-r", "--resolution", type=float, default=0.5, help="Grid resolution"
)
p.add_argument(
    "-d", "--dimension", type=float, default=23.5, help="Grid dimension (1D)"
)
p.add_argument(
    "-m", "--ligmap", type=str, default=f"{root}/files/ligmap", help="Ligand types file"
)
p.add_argument(
    "-o",
    "--output",
    type=str,
    default="grids",
    help="Output folder",
)
p.add_argument("--dx", action="store_true", help="Output grids as DX files")

args = p.parse_args()

# limbolgrid typer
typer = molgrid.FileMappedGninaTyper(args.ligmap)

mol = next(pybel.readfile("sdf", args.mol_sdf))
scaffold = next(pybel.readfile("sdf", args.scaffold_sdf))

gmol, cmol = mol_to_grid(mol, args.dimension, args.resolution, typer)
gscaffold, cscaffold = mol_to_grid(scaffold, args.dimension, args.resolution, typer, c=cmol)

gdiff = gmol - gscaffold

# !!! type_names argument works with the densitymatch container
# !!! 

print("Writing DX for molecule...")
molgrid.write_dx_grids(
    prefix=str(os.path.join(args.output, "mol")),
    type_names=typer.get_type_names(),
    grid=gmol[0].cpu(),
    center=cmol,
    resolution=args.resolution,
    scale=1.0,
)

print("Writing DX for scaffold...")
molgrid.write_dx_grids(
    prefix=str(os.path.join(args.output, "scaffold")),
    type_names=typer.get_type_names(),
    grid=gscaffold[0].cpu(),
    center=cmol,
    resolution=args.resolution,
    scale=1.0,
)

print("Writing DX for difference...")
molgrid.write_dx_grids(
    prefix=str(os.path.join(args.output, "diff")),
    type_names=typer.get_type_names(),
    grid=gdiff[0].cpu(),
    center=cmol,
    resolution=args.resolution,
    scale=1.0,
)