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
Calculate interactions using ODDT (https://oddt.readthedocs.io/).

Example usage:
    $ python -m pipelines.xchem.calc_interactions -p ../../data/mpro/Mpro-x0387_0.pdb -i ../../data/mpro/hits-17.sdf.gz -o output

The interactions that are searched for are:
* Hydrogen bond
* Halogen bond
* Hydrophobic
* Salt bridge
* Pi stacking
* Pi cation
See the ODDT docs for more details.

We introduce on additional concept, the 'canonical' interaction at the protein.
This allows the same interaction to be identified when comparing structures of the same protein.
The canonical representation comprise the protein residue sometimes followed in the case of H-bond interactions by a
definition of which part of that residue is involved (BN: backbone N, BO: backbone O, SC, sidechain). Other types of
interaction do not need this distinction so just use the residue type and number.

In the case of hydrophobic interactions there can be multiple interactions between a ligand and a residue. The information
recorded is only for the shortest of these interactions.

The different interaction types are written as data fields with a syntax like:

GLN189SC [10.922, 4.777, 26.595] [12.074, 5.388, 23.622] 3.246 3
HIS41SC [12.815, -4.64, 20.946] [13.769, -1.898, 21.404] 2.939 15

The tokens are:
1. The canonical site on the protein
2. The 3D location of the protein atom or pseudo atom as x, y, z coordinates inside square brackets
3. The 3D location of the ligand atom or pseudo atom as x, y, z coordinates inside square brackets
4. The distance between the atoms or pseudo atoms
5. Atom number (zero based) of the ligand atom. This is missing where the ligand part is a pseudo atom.

Normally the interaction is between 2 atoms, but in some cases like Pi interactions it is to an aromatic ring that is
represented as a pseudo atom at the centre of the ring.

You can also specify 'key interactions' that are counted. e.g in the case of the above H-bond interactions you can use
the `--key-hbond GLN189SC HIS41SC' option to count those specific interactions. One or more values are allowed. Similar
things for other types of interaction using the `--key-*` options.

Optionally, predicted binding scores using RDScore and NNScore can be generated. To calculate these you must specify the
appropriate pickle containing the model. Default models are provided based on PDB Bind data and can be specified as
follows:
* --rfscore RFScore_v1_pdbbind2016.pickle
* --nnscore NNScore_pdbbind2016.pickle

To regenerate these models use the build_oddt_models module.

Typically you specify a single protein and a set of ligands (e.g. as a SDF) to work on. Each ligand is calculated against
that protein. However, you can also specify multiple proteins as the --protein option allows multiple values.
In this case the proteins are used in turn (first ligand uses the first protein, second ligand uses the second protein ...)
until the proteins run out after which the last one continues to be used. In the usual case just specify a single protein
and it will be used for all ligands.

The following data fields will be generated. Not all will be present depending of the options used and whether particular
interactions are found:
* HydrophobicInteractions - details of the H-bond interactions
* HalogenBondInteractions - details of the halogen bond interactions
* HydrophobicInteractions - details of the hydrophobic interactions
* SaltBridgeInteractions - details of the salt bridge interactions
* PiStackingInteractions - details of the Pi stacking interactions
* PiCationInteractions - details of the Pi cation interactions
* NumTotalInteractions - total number of interactions found
* NumKeyInteractions - number of key interactions found
* KeyInteractions - the key interactions that were found
* vina_* - a whole series of AutoDock VINA scores that are generated when running NNScore
"""

from __future__ import print_function
import argparse, traceback
import json

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils
import oddt
from oddt import toolkit, spatial
from oddt.scoring import functions
import numpy as np

from . import interactions


# start function definitions #########################################

def generate_js(type, values):
    content = []
    for v in values:
        js = "  ['%s', '%s', [%s, %s, %s], [%s, %s, %s]]," % (
            type, v[0], v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2])
        content.append(js)
    return "\n".join(content)


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


def get_canonical_residue(atom):
    return atom['resname'] + str(atom['resnum'])


def determine_protein_format(protein_file, protein_format):
    if protein_format:
        return protein_format
    elif protein_file.endswith('.pdb') or protein_file.endswith('.pdb.gz'):
        return 'pdb'
    elif protein_file.endswith('.mol2') or protein_file.endswith('.mol2.gz'):
        return 'mol2'
    else:
        utils.log("Can't determine file format, so assuming pdb. Please use the --protein-format parameter")
        return 'pdb'


def read_next_protein(proteins, format, previous, index, keep_hs=False):
    if previous and index >= len(proteins):
        return previous
    protein = next(toolkit.readfile(format, proteins[index], removeHs=not keep_hs))
    if not protein:
        raise ValueError('Unable to read protein')
    else:
        utils.log('Read protein', index + 1)
        protein.protein = True
        protein.removeh()
        return protein


NNSCORES = []
RFSCORES = []


def init_nnscore(files):
    global NNSCORES
    for file in files:
        NNSCORES.append(functions.nnscore().load(file))


def init_rfscore(files):
    global RFSCORES
    for file in files:
        RFSCORES.append(functions.rfscore().load(file))


def calc_nnscore(protein, ligand):
    for model in NNSCORES:
        model.set_protein(protein)
        model.predict_ligand(ligand)


def calc_rfscore(protein, ligand):
    for model in RFSCORES:
        model.set_protein(protein)
        model.predict_ligand(ligand)


def process(protein_files, ligands, writer, key_inters, protein_format=None, filter_strict=False,
            exact_protein=False, exact_ligand=False, keep_hs_protein=False, keep_hs_ligand=False, report_file=None,
            compare_file=None, nnscores=None, rfscores=None, plecscores=None):
    pformat = determine_protein_format(protein_files[0], protein_format)
    utils.log('Protein format:', pformat)
    utils.log(len(protein_files), 'proteins specified')
    ligands = toolkit.readfile('sdf', ligands, removeHs=not keep_hs_ligand)

    if report_file:
        report_data = []
    else:
        report_data = None

    if compare_file:
        with open(compare_file, "r") as f:
            txt = f.read()
            compare_data = interactions.from_json(txt)
    else:
        compare_data = None

    if nnscores:
        print('Initialising NNSCORE')
        init_nnscore(nnscores)
    if rfscores:
        print('Initialising RFSCORE')
        init_rfscore(rfscores)

    total = 0
    count = 0
    errors = 0
    protein = None
    for ligand in ligands:
        # print('Processing ligand', total + 1)
        try:
            protein = read_next_protein(protein_files, pformat, protein, total, keep_hs=keep_hs_protein)
            if nnscores:
                calc_nnscore(protein, ligand)
            if rfscores:
                calc_rfscore(protein, ligand)
            inter_data = process_mol(protein, ligand, key_inters, total, filter_strict=filter_strict,
                                     exact_protein=exact_protein, exact_ligand=exact_ligand, compare_data=compare_data)
            if report_data is not None:
                report_data.append(inter_data)

            # write the RDKit mol
            writer.write(ligand.Mol)
            count += 1
        except:
            errors += 1
            traceback.print_exc()
        finally:
            total += 1

    # print(json.dumps(report_data, cls=interactions.InteractionEncoder))

    if report_data:
        with open(report_file, 'w') as report:
            json.dump(report_data, report, cls=interactions.InteractionEncoder)

    return count, errors


def process_mol(protein, mol, key_inters, index, filter_strict=False,
                exact_protein=False, exact_ligand=False, compare_data=None):
    mol_key_inters = []
    inter_data = {}
    num_inters = 0
    all_scores = {}

    inter_set = interactions.InteractionSet(index, None, [])

    inter_data['Index'] = index
    if '_Name' in mol.data:
        molname = mol.data['_Name']
        inter_data['LigandName'] = molname
        inter_set.ligand_name = molname

    # handle H-bond interactions
    inters = calc_hydrogen_bond_interactions(protein, mol, key_inters.get(interactions.I_TYPE_HBOND, None),
                                             mol_key_inters, filter_strict=filter_strict,
                                             exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_HBOND] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_HBOND] = data
        inter_set.add(inters)

    # handle hydrophobic interactions
    inters = calc_hydrophobic_interactions(protein, mol, key_inters.get(interactions.I_TYPE_HYDROPHOBIC, None),
                                           mol_key_inters)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_HYDROPHOBIC] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_HYDROPHOBIC] = data
        inter_set.add(inters)

    # handle salt bridge interactions
    inters = calc_salt_bridge_interactions(protein, mol, key_inters.get(interactions.I_TYPE_SALT_BRIDGE, None),
                                           mol_key_inters, exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_SALT_BRIDGE] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_HYDROPHOBIC] = data
        inter_set.add(inters)

    # handle pi stacking interactions
    inters = calc_pi_stacking_interactions(protein, mol, key_inters.get(interactions.I_TYPE_PI_STACKING, None),
                                           mol_key_inters, filter_strict=filter_strict)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_PI_STACKING] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_PI_STACKING] = data
        inter_set.add(inters)

    # handle pi cation interactions
    inters = calc_pi_cation_interactions(protein, mol, key_inters.get(interactions.I_TYPE_PI_CATION, None),
                                         mol_key_inters, filter_strict=filter_strict,
                                         exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_PI_CATION] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_PI_CATION] = data
        inter_set.add(inters)

    # handle halogen bond interactions
    inters = calc_halogen_bond_interactions(protein, mol, key_inters.get(interactions.I_TYPE_HALOGEN, None),
                                            mol_key_inters, filter_strict=filter_strict)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[interactions.I_NAME_HALOGEN] = data
        num_inters += len(inters.interactions)
        inter_data[interactions.I_NAME_HALOGEN] = data
        inter_set.add(inters)

    mol.data['NumTotalInteractions'] = str(num_inters)
    mol.data['NumKeyInteractions'] = str(len(mol_key_inters))
    if mol_key_inters:
        keydata = []
        for key_inter in mol_key_inters:
            score = all_scores.get(key_inter, 0)
            mol.data[key_inter] = score
            keydata.append(key_inter + ' ' + str(score))
        mol.data['KeyInteractions'] = '\n'.join(keydata)

    return inter_set


def compare_interactions(inter_type, compare_to, all_scores):
    if inter_type and compare_to:
        for inter in inter_type.interactions:
            # print('  looking at', inter_type.interaction_type, inter.canon_residue, inter.ligand_pos)
            for compare_set in compare_to:
                for compare_type in compare_set.interaction_types:
                    if compare_type.interaction_type == inter_type.interaction_type:
                        # print('    analysing', compare_set.ligand_name, compare_type.interaction_type, len(compare_type.interactions))
                        for compare_inter in compare_type.interactions:
                            # print('      checking', inter.canon_residue, compare_inter.canon_residue)
                            if inter.canon_residue == compare_inter.canon_residue:
                                # print('     ', compare_inter.canon_residue, compare_set.ligand_name)
                                sc = inter.compare_and_store(compare_inter, compare_set.ligand_name)

            if inter.score_value:
                k = inter_type.interaction_type[:-11] + ':' + inter.canon_residue
                # print('Found key interaction', k, inter.score_value)
                if k in all_scores:
                    # keep the best score
                    if inter.score_value > all_scores[k]:
                        all_scores[k] = inter.score_value
                else:
                    all_scores[k] = inter.score_value


def calc_hydrogen_bond_interactions(protein, mol, key_inters_defs, mol_key_inters,
                                    filter_strict=False, exact_protein=False, exact_ligand=False):
    """ Calculate H-bond interactions

    Parameters:
    protein (Molecule): The protein
    mol (Molecule): The ligand to test
    key_inters_defs (List): Optional list of key H-bond interactions
    mol_key_inters (List): List to add the key H-bond interactions to
    filter_strict (Bool): Whether to use strict matching

    Returns:
    InteractionType: The interactions
    """

    # first apply a fix that's needed to handle the mis-assignment of donor atoms for a particular tautomeric
    # form: nc[n;H1]. See https://groups.google.com/g/oddt/c/fqzmhSqprw8/m/nmaaUlCDAgAJ
    # mol.atom_dict.setflags(write=True)
    # matches = oddt.toolkit.Smarts('[n;H0]').findall(mol)
    # for (idx,) in matches:
    #     # idx assumes 0 based indexing e.g. RDKit. OBabel uses 1 based indexing.
    #     mol.atom_dict['isdonor'][idx] = False
    # mol.atom_dict.setflags(write=False)
    # end of fix

    protein_atoms, ligand_atoms, strict = oddt.interactions.hbonds(protein, mol,
                                                                   mol1_exact=exact_protein, mol2_exact=exact_ligand)
    inters = {}

    for p, l, s in zip(protein_atoms, ligand_atoms, strict):
        if s or not filter_strict:
            c = get_canonical_hbond(p)
            dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
            p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
            l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
            t = interactions.Interaction(c, p_coords, l_coords, dist, l['id'].item())
            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(interactions.I_TYPE_HBOND + ':' + c)
    if inters:
        # print('  found', len(inters), 'h-bonds')
        return interactions.InteractionType(interactions.I_NAME_HBOND, list(inters.values()))
    else:
        return None


def calc_hydrophobic_interactions(protein, mol, key_inters_defs, mol_key_inters):
    """ Calculate hydrophobic interactions.
    ODDT generates hydrophobic interactions for all pairs of atoms that are within the specification which results in
    multiple interactions between the same hydrophobic region of the ligand and protein.
    We reduce those down to a single interaction, the one that is the shortest.
    ODDT also seems to generate hydrophobic interactions that are complete duplicates, and the de-duplication process
    removes these as well.

    Parameters:
    protein (Molecule): The protein
    mol (Molecule): The ligand to test
    key_inters_defs (List): Optional list of key H-bond interactions
    mol_key_inters (List): List to add the key H-bond interactions to

    Returns:
    InteractionType: The interactions
    """
    inters = {}
    protein_atoms, ligand_atoms = oddt.interactions.hydrophobic_contacts(protein, mol)
    for p, l in zip(protein_atoms, ligand_atoms):
        c = get_canonical_residue(p)
        dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
        p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
        l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
        t = interactions.Interaction(c, p_coords, l_coords, dist, l['id'].item())
        if c in inters:
            current_value = inters[c]
            if dist < current_value.distance:
                inters[c] = t
        else:
            inters[c] = t
            if key_inters_defs and c in key_inters_defs:
                mol_key_inters.append(interactions.I_TYPE_HYDROPHOBIC + ':' + c)

    if inters:
        # print('  found', len(inters), 'hydrophobics')
        return interactions.InteractionType(interactions.I_NAME_HYDROPHOBIC, list(inters.values()))
    else:
        return None


def calc_salt_bridge_interactions(protein, mol, key_inters_defs, mol_key_inters, exact_protein=False, exact_ligand=False):
    inters = {}
    protein_atoms, ligand_atoms = oddt.interactions.salt_bridges(protein, mol, mol1_exact=exact_protein, mol2_exact=exact_ligand)
    for p, l in zip(protein_atoms, ligand_atoms):
        c = get_canonical_residue(p)
        dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
        p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
        l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
        t = interactions.Interaction(c, p_coords, l_coords, dist, l['id'].item())
        if c in inters:
            current_value = inters[c]
            if dist < current_value.distance:
                inters[c] = t
        else:
            inters[c] = t
            if key_inters_defs and c in key_inters_defs:
                mol_key_inters.append(interactions.I_TYPE_SALT_BRIDGE + ':' + c)

    if inters:
        return interactions.InteractionType(interactions.I_NAME_SALT_BRIDGE, list(inters.values()))
    else:
        return None


def calc_pi_stacking_interactions(protein, mol, key_inters_defs, mol_key_inters, filter_strict=False):
    protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular = oddt.interactions.pi_stacking(protein, mol)
    inters = {}
    for p, l, s1, s2 in zip(protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular):
        if (s1 or s2) or not filter_strict:
            c = get_canonical_residue(p)
            dist = spatial.distance(np.array([p['centroid']]), np.array([l['centroid']]))[0][0]
            p_coords = (p['centroid'][0].item(), p['centroid'][1].item(), p['centroid'][2].item())
            l_coords = (l['centroid'][0].item(), l['centroid'][1].item(), l['centroid'][2].item())
            t = interactions.Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(interactions.I_TYPE_PI_STACKING + ':' + c)

    if inters:
        return interactions.InteractionType(interactions.I_NAME_PI_STACKING, list(inters.values()))
    else:
        return None


def calc_pi_cation_interactions(protein, mol, key_inters_defs, mol_key_inters,
                                filter_strict=False, exact_protein=False, exact_ligand=False):
    """
    Pi-cation calculations are directional so are run in both directions (protein->ligand and ligand-> protein).
    Hence this method runs the pi_cation() calculation twice.
    :param protein:
    :param mol:
    :param key_inters_defs:
    :param mol_key_inters:
    :param filter_strict:
    :param exact_protein:
    :param exact_ligand:
    :return:
    """
    inters = {}
    # first treat the protein as the pi and the ligand as the cation
    rings, cation, strict = oddt.interactions.pi_cation(protein, mol, cation_exact=exact_ligand)
    for ring, cat, s in zip(rings, cation, strict):
        if s or not filter_strict:
            dist = spatial.distance(np.array([ring['centroid']]), np.array([cat['coords']]))[0][0]
            c = get_canonical_residue(ring)
            p_coords = (ring['centroid'][0].item(), ring['centroid'][1].item(), ring['centroid'][2].item())
            l_coords = (cat['coords'][0].item(), cat['coords'][1].item(), cat['coords'][2].item())
            t = interactions.Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(interactions.I_TYPE_PI_CATION + ':' + c)

    # now with the ligand as the pi and the protein as the cation
    rings, cation, strict = oddt.interactions.pi_cation(mol, protein, cation_exact=exact_protein)
    for ring, cat, s in zip(rings, cation, strict):
        if s or not filter_strict:
            dist = spatial.distance(np.array([ring['centroid']]), np.array([cat['coords']]))[0][0]
            c = get_canonical_residue(ring)
            p_coords = (cat['coords'][0].item(), cat['coords'][1].item(), cat['coords'][2].item())
            l_coords = (ring['centroid'][0].item(), ring['centroid'][1].item(), ring['centroid'][2].item())
            t = interactions.Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(interactions.I_TYPE_PI_CATION + ':' + c)


    if inters:
        return interactions.InteractionType(interactions.I_NAME_PI_CATION, list(inters.values()))
    else:
        return None


def calc_halogen_bond_interactions(protein, mol, key_inters_defs, mol_key_inters, filter_strict=False):
    protein_atoms, ligand_atoms, strict = oddt.interactions.halogenbonds(protein, mol)
    inters = {}
    for p, l, s in zip(protein_atoms, ligand_atoms, strict):
        if s or not filter_strict:
            c = get_canonical_residue(p)
            dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
            p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
            l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
            t = interactions.Interaction(c, p_coords, l_coords, dist, l['id'].item())
            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(interactions.I_TYPE_HALOGEN + ':' + c)

    if inters:
        return interactions.InteractionType(interactions.I_NAME_HALOGEN, list(inters.values()))
    else:
        return None


### start main execution #########################################

def main():
    # Example usage
    # python -m pipelines.xchem.calc_interactions -p ../../data/mpro/Mpro-x0387_0.pdb -i ../../data/mpro/hits-17.sdf.gz -o output

    parser = argparse.ArgumentParser(description='Calculate interactions')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-p', '--protein', nargs='*', help="File with protein (PDB or MOL2 format")
    # NOTE reading mol2 format seems to be problematical.
    parser.add_argument('-pf', '--protein-format', choices=['pdb', 'mol2'], help="Protein file format")
    parser.add_argument('--strict', action='store_true', help='Strict filtering')
    parser.add_argument('--exact-protein', action='store_true', help='Exact matching of hydrogens and charges for protein')
    parser.add_argument('--exact-ligand', action='store_true', help='Exact matching of hydrogens and charges for ligand')
    parser.add_argument('--keep-hs-protein', action='store_true', help='Keep hydrogens on the protein')
    parser.add_argument('--keep-hs-ligand', action='store_true', help='Keep hydrogens on the ligand')
    parser.add_argument('--key-hbond', nargs='*', help='List of canonical H-bond interactions to count')
    parser.add_argument('--key-hydrophobic', nargs='*', help='List of canonical hydrophobic interactions to count')
    parser.add_argument('--key-salt-bridge', nargs='*', help='List of canonical salt bridge interactions to count')
    parser.add_argument('--key-pi-stacking', nargs='*', help='List of canonical pi stacking interactions to count')
    parser.add_argument('--key-pi-cation', nargs='*', help='List of canonical pi cation interactions to count')
    parser.add_argument('--key-halogen', nargs='*', help='List of canonical halogen bond interactions to count')
    parser.add_argument('--rfscores', nargs='*', help="Pickle(s) for RFScore model e.g. RFScore_v1_pdbbind2016.pickle")
    parser.add_argument('--nnscores', nargs='*', help="Pickle(s) for NNScore model e.g. NNScore_pdbbind2016.pickle")
    parser.add_argument('--plecscores', nargs='*',
                        help="Pickle(s) for PLECScore model e.g. PLEClinear_p5_l1_pdbbind2016_s65536.pickle")
    parser.add_argument('-r', '--report-file', help="File for the report")
    parser.add_argument('-c', '--compare', help="Compare interactions with this report")

    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')

    args = parser.parse_args()
    utils.log("Calculate interactions Args: ", args)

    key_inters = {}
    if args.key_hbond:
        key_inters[interactions.I_TYPE_HBOND] = args.key_hbond
    if args.key_hydrophobic:
        key_inters[interactions.I_TYPE_HYDROPHOBIC] = args.key_hydrophobic
    if args.key_salt_bridge:
        key_inters[interactions.I_TYPE_SALT_BRIDGE] = args.key_salt_bridge
    if args.key_pi_stacking:
        key_inters[interactions.I_TYPE_PI_STACKING] = args.key_pi_stacking
    if args.key_pi_cation:
        key_inters[interactions.I_TYPE_PI_CATION] = args.key_pi_cation
    if args.key_halogen:
        key_inters[interactions.I_TYPE_HALOGEN] = args.key_halogen

    source = "calc_interactions.py using ODDT"
    datasetMetaProps = {"source": source, "description": "Calculate interactions using ODDT"}

    clsMappings = {}
    fieldMetaProps = []
    clsMappings[interactions.I_NAME_HBOND] = "java.lang.String"
    clsMappings[interactions.I_NAME_HALOGEN] = "java.lang.String"
    clsMappings[interactions.I_NAME_HYDROPHOBIC] = "java.lang.String"
    clsMappings[interactions.I_NAME_SALT_BRIDGE] = "java.lang.String"
    clsMappings[interactions.I_NAME_PI_STACKING] = "java.lang.String"
    clsMappings[interactions.I_NAME_PI_CATION] = "java.lang.String"
    clsMappings['NumTotalInteractions'] = "java.lang.Integer"
    clsMappings['NumKeyInteractions'] = "java.lang.Integer"
    clsMappings['KeyInteractions'] = "java.lang.String"
    fieldMetaProps.append({
        "fieldName": interactions.I_NAME_HBOND,
        "values": {"source": source, "description": "Hydrogen bond interactions"},
        "fieldName": interactions.I_NAME_HALOGEN,
        "values": {"source": source, "description": "Halogen bond interactions"},
        "fieldName": interactions.I_NAME_HYDROPHOBIC,
        "values": {"source": source, "description": "Hydrophobic interactions"},
        "fieldName": interactions.I_NAME_SALT_BRIDGE,
        "values": {"source": source, "description": "Salt bridge interactions"},
        "fieldName": interactions.I_NAME_PI_STACKING,
        "values": {"source": source, "description": "Pi stacking interactions"},
        "fieldName": interactions.I_NAME_PI_CATION,
        "values": {"source": source, "description": "Pi cation interactions"}
    })

    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                                                                  'calc_interactions', args.outformat,
                                                                  valueClassMappings=clsMappings,
                                                                  datasetMetaProps=datasetMetaProps,
                                                                  fieldMetaProps=fieldMetaProps,
                                                                  compress=not args.no_gzip)

    # this does the processing
    count, errors = process(args.protein, args.input, writer, key_inters,
                            protein_format=args.protein_format, filter_strict=args.strict,
                            exact_protein=args.exact_protein, exact_ligand=args.exact_ligand,
                            keep_hs_protein=args.keep_hs_protein, keep_hs_ligand=args.keep_hs_ligand,
                            report_file=args.report_file, compare_file=args.compare,
                            rfscores=args.rfscores, nnscores=args.nnscores, plecscores=args.plecscores)
    utils.log('Processing complete.', count, 'records processed.', errors, 'errors')

    inputs_file.close()
    writer.flush()
    writer.close()
    output.close()
    #
    if args.metrics:
        utils.write_metrics(output_base, {'__InputCount__': total, '__OutputCount__': count, '__ErrorCount__': errors,
                                          'ODDTInteraction': count})


if __name__ == "__main__":
    main()
