from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MCS import FindMCS
import os
import utils
import argparse


def guess_substruct(mol_one, mol_two):
    """Code to find the substructure between two molecules."""
    return Chem.MolToSmiles(Chem.MolFromSmarts(FindMCS([mol_one,mol_two],completeRingsOnly=True,matchValences=True).smarts))

def generate_conformers(my_mol, NumOfConf, ref_mol, outputfile, coreSubstruct):
    # Find the MCS if not given
    if not coreSubstruct:
        coreSubstruct = guess_substruct(my_mol,ref_mol)

    # Creating core of reference ligand #
    core1 = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(ref_mol,
                                                               Chem.MolFromSmiles(coreSubstruct)),
                                     Chem.MolFromSmiles('*'))
    core1.UpdatePropertyCache()

    # Generate conformers with constrained embed
    conf_lst = []
    for i in (xrange(NumOfConf)):
        conf_lst.append(Chem.AddHs(my_mol))
        AllChem.ConstrainedEmbed(conf_lst[i], core1, randomseed=i)
        outputfile.write(conf_lst[i])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RDKit constrained conformer generator')
    parser.add_argument('-n', '--num', type=int, default=10, help='number of conformers to generate')
    parser.add_argument('-r', '--refmol', help="Reference molecule file")
    parser.add_argument('-c', '--core_smi', help='Core substructure. If not specified - guessed using MCS', default='')
    utils.add_default_io_args(parser)


    args = parser.parse_args()
    # Get the reference molecule
    ref_mol_input, ref_mol_suppl = utils.default_open_input(args.refmol, args.refmol)
    counter = 0
    # Takes the last one in the SUPPL - sho
    for mol in ref_mol_suppl:
        ref_mol = mol
        counter+=1
    if counter >1:
        raise ValueError("Only one molecule should be given as reference. "+str(counter)+" were given.")

    # Get the molecules
    input, suppl = utils.default_open_input(args.input, args.informat)
    output, writer, output_base = utils.default_open_output(args.output, "rxn_maker", args.outformat)

    for mol in suppl:
        generate_conformers(mol, args.num, ref_mol, writer, args.core_smi)