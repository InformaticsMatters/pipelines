from rdkit import Chem

def load_data():
    """
    Function to load in the data for a generic cell
    :return: The data as an object - that can be used in the Ipython Notebook session.
    """
    # Input ones
    mol_input_one = "/data/input_mols.sdf"
    # For example we define other inputs - and these are just transferred here
    mol_input_two = "/data/input_mols_two.sdf"
    return Chem.SDMolSupplier(mol_input_one)