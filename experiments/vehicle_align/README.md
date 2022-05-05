# VEHICLe

## Notes

The installation of `pyvista` and `ipyvtklink` using `mamba` downgrades the `pymol-open-source` installation from `2.5` to `2.4`.

Using version `2.4` of `pymol-open-source` causes the `render.py` script to fail with the following error:

```text
 PyMOL not running, entering library mode (experimental)
terminate called after throwing an instance of 'pymol::ill_informed_image'
  what():  Image Construction ill-informed.
[1]    15257 abort (core dumped)  python render.py
```

Using version `2.5` of `pymol-open-source` produces the rendered images without problems.
