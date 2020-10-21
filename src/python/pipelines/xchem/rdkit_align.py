from rdkit import Chem
from rdkit.Chem import AllChem

def find_calphas(protein):
    calphas = {}
    for atom in protein.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        moninfo = atom.GetMonomerInfo()
        resnum = resinfo.GetResidueNumber()
        atomname = moninfo.GetName().strip()
        if 'CA' == atomname:
            calphas[resnum] = atom.GetIdx()
    return calphas

def align_calphas(probe, reference):

    ref_calphas = find_calphas(reference)
    print('Found', len(ref_calphas), 'CAs')
    prb_calphas = find_calphas(probe)
    print('Found', len(prb_calphas), 'CAs')
    atom_map = []
    for resnum, idx in prb_calphas.items():
        if resnum in ref_calphas:
            atom_map.append((idx, ref_calphas[resnum]))
        else:
            print('WARNING: residue', resnum, 'not found in reference')

    print('Mapped', len(atom_map), 'atoms')
    rmsd = AllChem.AlignMol(probe, reference, atomMap=atom_map)

    print('RMSD:', rmsd)

def extract_ligand(protein, resname):
    mol = Chem.RWMol(protein)
    atoms_to_delete = []
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        if resinfo.GetResidueName().strip() != resname:
            atoms_to_delete.append(atom.GetIdx())
    print('Deleting', len(atoms_to_delete), 'atoms')
    for idx in reversed(atoms_to_delete):
        mol.RemoveAtom(idx)
    return mol

def main():
    reference = Chem.MolFromPDBFile('hits23_complex_init_0.pdb')
    probe = Chem.MolFromPDBFile('hits23_complex_mini_0.pdb')
    align_calphas(probe, reference)
    Chem.MolToPDBFile(probe, 'hits23_complex_algn_0.pdb')

if __name__ == "__main__":
    main()
