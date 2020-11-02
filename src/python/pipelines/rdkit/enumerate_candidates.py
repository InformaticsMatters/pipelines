#!/usr/bin/env python

# Copyright 2020 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Reads a file of SMILES and enumerates undefined chiral centres, tautomers and charge states.
Writes results as SMILES.
"""
import sys, argparse, traceback

from pipelines_utils import utils

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField


from molvs.tautomer import TautomerEnumerator

from . import prepare_3d
from pipelines.dimorphite.dimorphite_dl import run_with_mol_list

def gen_tautomers(mols, enumerator):
    tautomers = []
    for m in mols:
        ts = enumerator.enumerate(m)
        tautomers.extend(ts)
    return tautomers

def gen_charges(mols, min_ph, max_ph, min_charge, max_charge,  num_charges):
    protonated_mols = run_with_mol_list(mols, min_ph=min_ph, max_ph=max_ph)
    wanted_mols = []
    for mol in protonated_mols:
        if num_charges is not None:
            nc = 0
            for atom in mol.GetAtoms():
                if atom.GetFormalCharge() > 0:
                    nc += 1
            if nc > num_charges:
                continue

        if min_charge is not None or max_charge is not None:
            charge = Chem.GetFormalCharge(mol)
            if min_charge is not None and charge < min_charge:
                continue
            if max_charge is not None and charge > max_charge:
                continue
        wanted_mols.append(mol)
    return wanted_mols

def add_molecules(mydict, mols):
    for mol in mols:
        add_molecule(mydict, mol)


def add_molecule(mydict, mol):
    smiles = Chem.MolToSmiles(mol)
    if smiles not in mydict:
        mydict[smiles] = mol

def execute(suppl, writer, enumerate_chirals=False,
            enumerate_charges=False, enumerate_tautomers=False,
            combinatorial=False,
            min_charge=None, max_charge=None, num_charges=None,
            min_ph=5.0, max_ph=9.0,
            gen3d=False, add_hydrogens=False, minimize=0, smiles_field=None):

    utils.log('Executing ...')

    if enumerate_tautomers or combinatorial:
        enumerator = TautomerEnumerator()

    if minimize:
        getFF = UFFGetMoleculeForceField

    count = 0
    total = 0
    errors = 0
    for mol in suppl:
        count += 1
        name = mol.GetProp('_Name')
        #utils.log('Processing mol', count)

        try:

            # use a dict as we want to keep the order
            enumerated_mols = {}
            add_molecule(enumerated_mols, mol)

            if enumerate_chirals:
                mols = prepare_3d.enumerate_undefined_chirals(mol)
                add_molecules(enumerated_mols, mols)

            if combinatorial:
                tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                protonated_mols = gen_charges(tautomers, min_ph, max_ph, min_charge, max_charge, num_charges)
                add_molecules(enumerated_mols, protonated_mols)
            else:
                if enumerate_tautomers:
                    tautomers = gen_tautomers(enumerated_mols.values(), enumerator)
                if enumerate_charges:
                    protonated_mols = gen_charges(enumerated_mols.values(), min_ph, max_ph, min_charge, max_charge, num_charges)
                if enumerate_tautomers:
                    add_molecules(enumerated_mols, tautomers)
                if enumerate_charges:
                    add_molecules(enumerated_mols, protonated_mols)


            total += len(enumerated_mols)

            for m in enumerated_mols.values():

                if smiles_field:
                    m.SetProp(smiles_field, Chem.MolToSmiles(m))

                if gen3d:
                    if add_hydrogens:
                        m = Chem.AddHs(m)
                    m = prepare_3d.convert_to_3d(m)
                    if not m:
                        errors += 1
                        continue
                    if minimize:
                        prepare_3d.minimize_mol(m, getFF, minimize)

                m.SetProp('_Name', name)
                writer.write(m)

        except:
            errors += 1
            traceback.print_exc()


    return count, total, errors


### start main execution #########################################

def main():

    # Example:
    #   python -m pipelines.rdkit.enumerate_candidates -i ../../data/nci10.smiles -o foo.smi --enumerate-charges --enumerate-charges --enumerate-tautomers

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate candidates')
    parser.add_argument('-i', '--input', help="Input file as SMILES")
    parser.add_argument('-o', '--output', help="Output file as SMILES or SDF")
    parser.add_argument('-d', '--delimiter', default='\t', help="Delimiter")
    parser.add_argument("--name-column", type=int, default=1, help="Column for name field (zero based)")
    parser.add_argument('-t', '--title-line', default=False, help="Do the files have title lines")
    parser.add_argument('--enumerate-charges', help='Enumerate charge forms', action='store_true')
    parser.add_argument('--enumerate-chirals', help='Enumerate undefined chiral centers', action='store_true')
    parser.add_argument('--enumerate-tautomers', help='Enumerate undefined chiral centers', action='store_true')
    parser.add_argument('--combinatorial', help='Combinatorial enumeration of charge and tautomers', action='store_true')
    parser.add_argument("--min-ph", type=float, default=5, help="Minimum pH for charge enumeration")
    parser.add_argument("--max-ph", type=float, default=9, help="Maximum pH for charge enumeration")
    parser.add_argument("--min-charge", type=int, help="Minimum charge of molecule to process")
    parser.add_argument("--max-charge", type=int, help="Maximum charge of molecule to process")
    parser.add_argument("--num-charges", type=int, help="Maximum number of atoms with a charge")
    parser.add_argument('--add-hydrogens', action='store_true', help='Include hydrogens in the output')

    parser.add_argument('--gen3d', action='store_true', help='Generate 3D coordinates (output as SDF)')
    parser.add_argument("--minimize", type=int, default=0, help="number of minimisation cycles when generating 3D molecules")
    parser.add_argument('--smiles-field', help='Add the SMILES as a field of this name (SDF output only)')

    args = parser.parse_args()
    utils.log("enumerate_candidates: ", args)

    # save the arguments
    input = args.input
    output = args.output
    delimiter = args.delimiter
    title_line = args.title_line
    name_column = args.name_column
    enumerate_charges = args.enumerate_charges
    enumerate_chirals = args.enumerate_chirals
    enumerate_tautomers = args.enumerate_tautomers
    combinatorial = args.combinatorial
    min_ph = args.min_ph
    max_ph = args.max_ph
    min_charge = args.min_charge
    max_charge = args.max_charge
    num_charges = args.num_charges
    add_hydrogens = args.add_hydrogens
    gen3d = args.gen3d
    minimize = args.minimize
    smiles_field = args.smiles_field


    # Dimporphite needs to use argparse with its own arguments, not messed up with our arguments
    # so we store the original args
    orig_sys_argv = sys.argv[:]

    # Remove all the parameters, keeping only the filename (first one) so that
    # dimorphite is unaware of any previous commandline parameters.
    sys.argv = sys.argv[:1]

    suppl = Chem.SmilesMolSupplier(input, delimiter=delimiter, titleLine=title_line, nameColumn=name_column)

    if gen3d:
        writer = Chem.SDWriter(output)
    else:
        writer = Chem.SmilesWriter(output, delimiter=delimiter, includeHeader=title_line)

    count, total, errors = execute(suppl, writer, enumerate_chirals=enumerate_chirals,
                                   enumerate_charges=enumerate_charges, enumerate_tautomers=enumerate_tautomers,
                                   combinatorial=combinatorial,
                                   min_charge=min_charge, max_charge=max_charge, num_charges=num_charges,
                                   min_ph=min_ph, max_ph=max_ph, add_hydrogens=add_hydrogens,
                                   gen3d=gen3d, minimize=minimize, smiles_field=smiles_field
                                   )

    utils.log(count, total, errors)

    writer.close()




if __name__ == "__main__":
    main()
