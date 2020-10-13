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
Prepare 3D molecules for docking.

This module can be used to enumerate undefined chiral centers and then generate a 3D representation
of the molecule. It is typically used to prepare ligands for docking.
There are also options for filtering out molecules with excess charges and various output options.
See the arguments in the main() method for full details,

Example usage:
    $ python -m pipelines.rdkit.prepare_3d -i ../../data/mpro/hits-17.sdf.gz -o output

"""

import argparse

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import *

def enumerateMol(mol, fragment):
    """
    Enumerate a single molecule
    :param mol:
    :param fragment The fragmentation method, 'hac' or 'mw'. If not specified the whole molecules is passed to Dimorphite
    :return:
    """

    if fragment:
        mol = mol_utils.fragment(mol, fragment)

    inputmol = []
    inputmol.append(mol)

    protonated_mols = run_with_mol_list(inputmol)
    return protonated_mols


def writeEnumeratedMols(src_mol, enum_mols, writer, index, name):
    """
    Write the enumerated molecule to the writer, incorporating the ID of the source molecule.
    :param src_mol:
    :param enum_mols:
    :param writer:
    :param index:
    :return:
    """
    errors = 0
    total = 0
    for mol in enum_mols:
        if mol is None:
            errors += 1
            continue
        add_src_mol_ref(src_mol, mol, index, name)
        writer.write(mol)
        total += 1

    return total, errors


def add_src_mol_ref(src_mol, target_mol, index, name):
    """
    Add the ID of the source molecule to the enumerated molecule as the field named SrcMolUUID.
    The ID is taken form the uuid field if it exists.
    The SrcMolIdx field is always set with the index of the source molecule in the input.
    The SrcMolName is taken from the _Name field if it exists.
    If the name parameter is specified it is added as the _Name field.
    :param src_mol:
    :param target_mol:
    :param index:
    :return:
    """

    parent = None
    if src_mol.HasProp('uuid'):
        parent = src_mol.GetProp('uuid')

    if parent:
        target_mol.SetProp('SrcMolUUID', parent)

    target_mol.SetIntProp('SrcMolIdx', index)
    if src_mol.HasProp('_Name'):
        n = src_mol.GetProp('_Name')
        if name:
            target_mol.SetProp('SrcMolName', n)
    if name:
        target_mol.SetProp('_Name', name)


def convert_to_3d(mol):
    molh = Chem.AddHs(mol, addCoords=True)
    ci = AllChem.EmbedMolecule(molh)
    if ci >= 0:
        return molh
    else:
        #print('No conformer')
        return None

def minimize_mol(mol, getForceField, cycles):
    ff = getForceField(mol, confId=-1)
    ff.Initialize()
    n = cycles
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        n -= 1

num_swap_success = 0

def enumerate_undefined_chirals(mol):
    """
    Enumerate undefined chiral centres that are points of substitution.
    Chiral centres that are already specifically defined are left untouched.
    :param mol: The mol to enumerate
    :return:
    """

    global num_swap_success

    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    #utils.log('  Chiral centers:', chiral_centers, 'Free atom counts:', free_counts)

    chiral_targets = []
    for i, cc in enumerate(chiral_centers):
        if cc[1] == '?':  # only enumerate if the stereo is undefined
            chiral_targets.append(cc[0])

    if chiral_targets:
        chiral_swaps = []
        mol_cw = Chem.RWMol(mol)
        mol_cw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        chiral_swaps.append(mol_cw)
        mol_ccw = Chem.RWMol(mol)
        mol_ccw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        chiral_swaps.append(mol_ccw)
        for j in range(1, len(chiral_targets)):
            new_enum_mols = []
            for m in chiral_swaps:
                nm1 = Chem.RWMol(m)
                nm1.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                new_enum_mols.append(nm1)
                nm2 = Chem.RWMol(m)
                nm2.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                new_enum_mols.append(nm2)
                chiral_swaps = new_enum_mols
                num_swap_success += len(chiral_swaps)
    else:
        chiral_swaps = [Chem.RWMol(mol)]
    #utils.log('  Mols to embed:', len(chiral_swaps))
    return chiral_swaps


def execute(molecules, writer, minimize=0, name_field=None,
            enumerate_chirals=0, smiles_field=None, include_hs=False, min_charge=None, max_charge=None, num_charges=None):
    count = 0
    total = 0
    errors = 0
    for mol in molecules:
        count += 1
        utils.log('Processing mol', count)
        if mol is None:
            continue
        if name_field:
            if mol.HasProp(name_field):
                name = mol.GetProp(name_field)
            else:
                name = ''
        else:
            name = Chem.MolToSmiles(mol)

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

        if enumerate_chirals:
            enum_mols = enumerate_undefined_chirals(mol)
            # utils.log('Enumerated', len(enum_mols), 'chirals')
        else:
            enum_mols = [mol]

        converted = []
        for m in enum_mols:
            # print('convert3d', Chem.MolToSmiles(m))
            mol3d = convert_to_3d(m)
            # print(mol3d)
            if not mol3d:
                utils.log('Failed to generate 3D for mol', count)
                errors += 1
                continue
            if minimize:
                minimize_mol(mol3d, UFFGetMoleculeForceField, minimize)
            converted.append(mol3d)

        # utils.log('Created', len(converted), '3D structures')
        enum_mols = converted

        if smiles_field:
            for m in enum_mols:
                m.SetProp(smiles_field, Chem.MolToSmiles(Chem.RemoveHs(m)))

        if not include_hs:
            adjusted_for_hs = []
            for m in enum_mols:
                adjusted_for_hs.append(Chem.RemoveHs(m))
            enum_mols = adjusted_for_hs

        utils.log('Created', len(enum_mols), '3D structures')
        t, e = writeEnumeratedMols(mol, enum_mols, writer, count, name)
        total += t
        errors += e
    return count, total, errors

### start main execution #########################################

def main():
    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate charges')
    parser.add_argument('-i', '--input',
                        help="Input file, if not defined the STDIN is used")
    parser.add_argument('-if', '--informat', choices=['sdf', 'json', 'smi'],
                        help="Input format. When using STDIN this must be specified.")
    parameter_utils.add_default_output_args(parser)
    parser.add_argument('--fragment-method', choices=['hac', 'mw'],
                        help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument("--minimize", type=int, default=0, help="number of minimisation cycles when generating 3D molecules")
    parser.add_argument('--enumerate-chirals', help='Enumerate undefined chiral centers', type=int, default=0)
    parser.add_argument('--smiles-field', help='Add the SMILES as a field of this name')
    parser.add_argument('--name-field', help='Use this field in the input as the name field in the output. If not specified the SMILES is used')
    parser.add_argument('--include-hydrogens', action='store_true', help='Include hydrogens in the output')
    parser.add_argument("--min-charge", type=int, help="Minimum charge of molecule to process")
    parser.add_argument("--max-charge", type=int, help="Maximum charge of molecule to process")
    parser.add_argument("--num-charges", type=int, help="Maximum number of atoms with a charge")

    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')

    args = parser.parse_args()
    utils.log("prepare_3d: ", args)

    # handle metadata
    source = "prepare_3d.py"
    datasetMetaProps = {"source": source, "description": "Enumerate undefined stereo isomers and generate 3D"}
    clsMappings = {
        "EnumChargesSrcMolUUID": "java.lang.String",
        "EnumChargesSrcMolIdx": "java.lang.Integer"
    }
    fieldMetaProps = [
        {"fieldName": "EnumChargesSrcMolUUID", "values": {"source": source, "description": "UUID of source molecule"}},
        {"fieldName": "EnumChargesSrcMolIdx", "values": {"source": source, "description": "Index of source molecule"}}
    ]

    oformat = utils.determine_output_format(args.outformat)
    if args.informat == 'smi':
        suppl = rdkit_utils.default_open_input_smiles(args.input, delimiter='\t', titleLine=False)
        input = None
    else:
        input, suppl = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                                                                 'enumerate_molecules', args.outformat,
                                                                 thinOutput=False, valueClassMappings=clsMappings,
                                                                 datasetMetaProps=datasetMetaProps,
                                                                 fieldMetaProps=fieldMetaProps,
                                                                 compress=not args.no_gzip)

    count, total, errors = execute(suppl, writer, minimize=args.minimize, enumerate_chirals=args.enumerate_chirals,
                                   smiles_field=args.smiles_field, include_hs=args.include_hydrogens,
                                   min_charge=args.min_charge, max_charge=args.max_charge, num_charges=args.num_charges)

    utils.log(count, total, errors)

    if input:
        input.close()
    writer.flush()
    writer.close()
    output.close()

    # re-write the metadata as we now know the size
    if oformat == 'json':
        utils.write_squonk_datasetmetadata(output_base, False, clsMappings, datasetMetaProps, fieldMetaProps, size=total)

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': count, '__OutputCount__': total, '__ErrorCount__': errors,
                                          'EnumerateChargesDimporphite': total})


if __name__ == "__main__":
    main()
