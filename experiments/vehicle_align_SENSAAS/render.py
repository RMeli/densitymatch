"""
Render top 3 alignment for all inhibitors for both CDK2 and BRD4
"""
from pymol import cmd
import os

def render(frag, lig):
    """
    Render alignment using PyMol API
    """
    cmd.reinitialize()
    cmd.load(frag, "fragmol")
    cmd.load(lig, "mol")
    cmd.orient("mol")
    cmd.center("mol")
    cmd.zoom("mol", 2.0)

    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_color", "black")
    cmd.set("ray_trace_mode", 3)
    cmd.set("antialias", 2)
    cmd.ray(1000, 1000)
    cmd.save(f"images/vehicle-{lig}.png")

    return cmd.ipython_image()


for system in ["BRD4"]:
    for lig in os.path.listdir(system):
        
        ffrag = os.path.join(system, lig, "bestsensaas.sdf")
        flig = os.path.join("../ligands/", system, f"{lig}.sdf")

        print(ffrag, flig)

        #render(ffrag, lig, i)
