from rdkit import Chem
def load_data():
    return Chem.SDMolSupplier("/data/input_mols.sdf")