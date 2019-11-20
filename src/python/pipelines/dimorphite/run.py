from rdkit import Chem
import dimorphite_dl

#suppl = [Chem.MolFromSmiles(s) for s in ["C[C@](F)(Br)CC(O)=O", "CCCCCN"]]
suppl = Chem.SDMolSupplier('Kinase_inhibs.sdf')
#suppl = Chem.SDMolSupplier('dhfr_standardized.sdf')

print(suppl)

protonated_mols = dimorphite_dl.run_with_mol_list(
    suppl,
    min_ph=5.0,
    max_ph=9.0,
)


print("Charged mols ------------------------------------------------")

for m in protonated_mols:
    if m:
        print(Chem.MolToSmiles(m) + " " + ",".join(m.GetPropNames()))
