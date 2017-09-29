from rdkit import Chem
def load_data():
    return Chem.SDMolsupplier("/data/input_mols.sdf")