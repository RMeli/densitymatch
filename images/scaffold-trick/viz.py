from pymol import cmd
import os
import glob


def viz(lig, dir="grids"):
    cmd.reinitialize("everything")
    
    sdfroot = "../../experiments/ligands/"
    ligpath = os.path.join(sdfroot, f"{lig}.sdf")
    scaffpath = os.path.join(sdfroot, f"{lig}_murcko.sdf")

    cmd.load(ligpath, "lig")
    cmd.load(scaffpath, "scaffold")

    cmd.remove("hydrogens")

    for f in glob.glob(os.path.join(dir, f"mol_*.dx")):
        namedx = f.replace("grids/", "")

        cmd.load(f, namedx)
        cmd.isomesh("m" + namedx, namedx, 0.5)

        if "Oxygen" in namedx:
            cmd.color("red", "m" + namedx)
        elif "Nitrogen" in namedx:
            cmd.color("blue", "m" + namedx)
        elif "Carbon" in namedx:
            if "Aromatic" in namedx:
                cmd.color("grey", "m" + namedx)
            else:
                cmd.color("white", "m" + namedx)
        else:
            cmd.color("purple", "m" + namedx)

    for f in glob.glob(os.path.join(dir, f"scaffold_*.dx")):
        namedx = f.replace("grids/", "")

        cmd.load(f, namedx)
        cmd.isomesh("m" + namedx, namedx, 0.5)

        if "Oxygen" in namedx:
            cmd.color("red", "m" + namedx)
        elif "Nitrogen" in namedx:
            cmd.color("blue", "m" + namedx)
        elif "Carbon" in namedx:
            if "Aromatic" in namedx:
                cmd.color("grey", "m" + namedx)
            else:
                cmd.color("white", "m" + namedx)
        else:
            cmd.color("purple", "m" + namedx)

    for f in glob.glob(os.path.join(dir, f"diff_*.dx")):
        namedx = f.replace("grids/", "")

        cmd.load(f, namedx)
        cmd.isomesh("m" + namedx, namedx, 0.5)

        if "Oxygen" in namedx:
            cmd.color("red", "m" + namedx)
        elif "Nitrogen" in namedx:
            cmd.color("blue", "m" + namedx)
        elif "Carbon" in namedx:
            if "Aromatic" in namedx:
                cmd.color("grey", "m" + namedx)
            else:
                cmd.color("white", "m" + namedx)
        else:
            cmd.color("purple", "m" + namedx)

cmd.extend("viz", viz)
