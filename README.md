# Density Match

Match [libmolgrid](https://github.com/gnina/libmolgrid)-generated densities using [sensaas](https://github.com/SENSAAS/sensaas).

## Pipeline

The goal is to align a molecule (or fragment) to a given `libmolgrid` density. This is achieved by transforming both the density and the given molecule into point clouds and aligning those.

### Fragment/Scaffold Placement and Molecular Growth

* [`libmolgrid`] Compute density for the molecule
* Compute point clouds from the densities (isosurfaces)
* [`SENSAAS`] Align point clouds
* Apply transformation to the molecule
* [`libmolgrid`] Re-compute density for the (aligned) molecule
* Compute difference between two densities
* [`liGAN`] Fit remaining density with atom
* [`liGAN`] Reconstruct full molecule

## Tools

### Convert SDF File to Point Cloud

A SDF file can be converted to a point cloud (`.pcd`, `.xyzrgb`, ...) using `molgrid_to_pcd.py`.
### View Point Cloud

One or multiple point clouds can be visualized with `view_pcd.py`.
### Convert Binary Point Cloud File to ASCII

A binary point cloud file can be converted into an ASCII file using `pcd_to_ascii.py`.
