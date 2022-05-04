"""
Render top 3 alignment for all inhibitors for both CDK2 and BRD4
"""

from rdkit import Chem
from pymol import cmd

import pickle
import pandas as pd

sorting_criteria = "gfit + hfit"


def render(df_mol, lig, idx=0):
    """
    Render alignment using PyMol API
    """
    mol = df_mol.iloc[idx].mol
    fragmol = df_mol.iloc[idx].fragmol

    with Chem.SDWriter("mol-tmp.sdf") as w:
        w.write(mol, confId=0)

    with Chem.SDWriter("fragmol-tmp.sdf") as w:
        w.write(fragmol, confId=1)

    cmd.reinitialize()
    cmd.load("fragmol-tmp.sdf", "fragmol")
    cmd.load("mol-tmp.sdf", "mol")
    cmd.orient("mol")
    cmd.center("mol")
    cmd.zoom("mol", 1.0)

    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_color", "black")
    cmd.set("ray_trace_mode", 3)
    cmd.set("antialias", 2)
    cmd.ray(1000, 1000)
    cmd.save(f"images/vehicle-{lig}-{idx}.png")

    return cmd.ipython_image()


for system in ["BRD4", "CDK2"]:

    with open(f"{system}-VEHICLe.pkl", "br") as fin:
        alignments = pickle.load(fin)

    # Convert nested index into tuples
    d = {}
    for outerKey, innerDict in alignments.items():
        for innerKey, values in innerDict.items():
            d[(outerKey, innerKey)] = values

    df = pd.DataFrame.from_dict(d)
    df = df.stack(level=0).swaplevel().sort_index()
    df.index.names = ["lig", "idx"]

    df["gfit + hfit"] = df["gfit"] + df["hfit"]

    for lig in list(df.index.levels[0]):
        df_mol = (
            df.query("lig == @lig")
            .sort_values(sorting_criteria, ascending=False)
            .reset_index()
        )

        for i in range(3):  # Top 3 alignment
            render(df_mol, lig, i)

# import matplotlib.image as mpimg
# import matplotlib.pyplot as plt

# fig, axs = plt.subplots(10, 3, figsize=(12, 10))

# for r, lig in enumerate(list(df.index.levels[0])):
#     for i in range(3):  # Top 3 alignment
#         render(df, lig, i)

#         print(r, i)

#         axs[r, i].imshow(mpimg.imread(f"images/vehicle-{lig}-{i}.png"))
#         axs[r, i].set_aspect("equal")
#         axs[r, i].set_xlabel(f"{lig.title()} (Rank {i})")
#         axs[r, i].set_ylabel("")
#         axs[r, i].set_xticks([])
#         axs[r, i].set_yticks([])


# plt.tight_layout()
# plt.savefig("vehicles-align-BRD4.png", dpi=600)
# plt.show()
