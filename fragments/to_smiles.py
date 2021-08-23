from rdkit.Chem import PandasTools

fname = "Enamine_DSI-Poised_Library.sdf"

df = PandasTools.LoadSDF(fname, smilesName="SMILES", molColName='mol')
df.drop(columns=["Catalog ID", "PlateID", "Well", "ID", "mol"], inplace=True)
df.rename(columns={"MW (desalted)": "MW"}, inplace=True)

df.to_csv("DSIpoised.csv", index=False)