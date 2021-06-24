# Density Match

Match [libmolgrid]()-generated densities using [sensaas]().

## Tools

### Convert SDF File to Point Cloud

A SDF file can be converted to a point cloud (`.pcd`, `.xyzrgb`, ...) using `molgrid_to_pcd.py`:
```text
usage: molgrid_to_pcd.py [-h] [-r RESOLUTION] [-d DIMENSION] [-m LIGMAP]
                         [-o OUTPUT] [--dx] [--ascii]
                         sdf

Visualize point cloud (PCD file) with Open3D

positional arguments:
  sdf                   SDF file

optional arguments:
  -h, --help            show this help message and exit
  -r RESOLUTION, --resolution RESOLUTION
                        Grid resolution
  -d DIMENSION, --dimension DIMENSION
                        Grid dimension (1D)
  -m LIGMAP, --ligmap LIGMAP
                        Ligand types file
  -o OUTPUT, --output OUTPUT
                        Output file (.pcd, .xyzrgb, ...)
  --dx                  Output grids as DX files
  --ascii               Output file in ASCII format
```

### View Point Cloud

One or multiple point clouds can be visualized with `view_pcd.py`:
```
usage: view_pcd.py [-h] files [files ...]

Visualize point cloud (PCD file) with Open3D

positional arguments:
  files       PCD files

optional arguments:
  -h, --help  show this help message and exit
```
### Convert Binary Point Cloud File to ASCII

A binary point cloud file can be converted into an ASCII file using `pcd_to_ascii.py`:
```text
usage: pcd_to_ascii.py [-h] [-o OUTPUT] file

Convert point cloud file to ASCII

positional arguments:
  file                  file to convert to ASCII

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output file
```