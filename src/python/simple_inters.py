import sys
import oddt
from oddt import interactions


if len(sys.argv) != 3:
    print("Usage: simple_inters.py protein.pdb ligand.mol")
    exit(1)

def get_canonical_hbond(atom):
    # print('classifying', atom['atomtype'], atom['isbackbone'], atom['isacceptor'], atom['isdonor'], atom['isdonorh'])
    res = atom['resname'] + str(atom['resnum'])
    if atom['isbackbone']:
        if atom['atomtype'] == 'N.am' or atom['atomtype'] == 'N.3':
            return res + 'BN'
        elif atom['atomtype'] == 'O.2':
            return res + 'BO'
        else:
            print('Unexpected H-bond atom', res, atom['atomtype'])
    else:
        return res + 'SC'

protein_pdbfile = sys.argv[1]
ligand_molfile = sys.argv[2]

exact_ligand = True

ligand = next(oddt.toolkit.readfile('sdf', ligand_molfile))
protein = next(oddt.toolkit.readfile('pdb', protein_pdbfile))
protein.protein = True

print('Protein:', protein_pdbfile)
print('Ligand:' + ligand_molfile)
print('Num protein/ligand atoms:', len(protein.atoms), len(ligand.atoms))
print('Exact ligand =', exact_ligand)

protein_atoms, ligand_atoms, strict = interactions.hbonds(protein, ligand, mol1_exact=False, mol2_exact=exact_ligand)
count = 0
for p, l, s in zip(protein_atoms, ligand_atoms, strict):
    count += 1
    print('  H-bond', get_canonical_hbond(p), '-', l['atomtype'], l['id'].item(), s)
print('Found', count, 'H-bond interactions')

protein_atoms, ligand_atoms = interactions.salt_bridges(protein, ligand, mol2_exact=exact_ligand)
count = 0
for p, l in zip(protein_atoms, ligand_atoms):
    count += 1
    print('  SaltBr', p['resname'] + str(p['resnum']), '-', l['atomtype'], l['id'].item())
print('Found', count, 'SaltBr interactions')

protein_atoms, ligand_atoms = oddt.interactions.hydrophobic_contacts(protein, ligand)
count = 0
for p, l in zip(protein_atoms, ligand_atoms):
    count += 1
    print('  Hphobe', p['resname'] + str(p['resnum']), '-', l['atomtype'], l['id'].item())
print('Found', count, 'Hphobe interactions')

protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular = oddt.interactions.pi_stacking(protein, ligand)
count = 0
for p, l, s1, s2 in zip(protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular):
    count += 1
    print('  PiStack', p['resname'] + str(p['resnum']), '-', s1, s2)
print('Found', count, 'pistack interactions')

count = 0
rings, cation, strict = oddt.interactions.pi_cation(protein, ligand, cation_exact=exact_ligand)
for ring, cat, s in zip(rings, cation, strict):
    count += 1
    print('  PiCation', ring['resname'] + str(ring['resnum']), 'protein-ligand -', s)
rings, cation, strict = oddt.interactions.pi_cation(ligand, protein, cation_exact=False)
for ring, cat, s in zip(rings, cation, strict):
    count += 1
    print('  PiCation', cat['resname'] + str(cat['resnum']), 'ligand-protein -', s)
print('Found', count, 'pication interactions')

protein_atoms, ligand_atoms, strict = oddt.interactions.halogenbonds(protein, ligand)
count = 0
for p, l, s in zip(protein_atoms, ligand_atoms, strict):
    count += 1
    print('  Halogen', p['resname'] + str(p['resnum']), '-', l['atomtype'], l['id'].item(), s)
print('Found', count, 'halogen interactions')
