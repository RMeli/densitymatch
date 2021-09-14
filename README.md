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

### Scripts

* [molgrid_to_pcd.py]: convert molecule into atomic density grid and subsequently convert the grid into a point cloud
* [score_pcd.py]: score alignment between point clouds and provide rototranslation for the optimal alignment
* [molgrid_diff.py]: compute difference between atomic density grids
* [fit_to_grid.py]: fit atoms to density grid (difference)

## Tools

### Convert SDF File to Point Cloud

A SDF file can be converted to a point cloud (`.pcd`, `.xyzrgb`, ...) using `molgrid_to_pcd.py`.
### View Point Cloud

One or multiple point clouds can be visualized with `view_pcd.py`.
### Convert Binary Point Cloud File to ASCII

A binary point cloud file can be converted into an ASCII file using `pcd_to_ascii.py`.


## Notes

### Running Experiments

Given [libmolgrid/#62](https://github.com/gnina/libmolgrid/issues/62) most scripts need to run within a Singularity container, where `libmolgrid` has been compiled from scratch. This extends to Jupyter Notebooks that require `libmolgrid`.

The `development/densitymatch.def` container allows to build a minimal Singularity container. The reconstruction (fitting) of atoms from atomic densities requires [liGAN](https://github.com/mattragoza/liGAN), which can be run within the full Singularity container built from `development/ligan.def`.
