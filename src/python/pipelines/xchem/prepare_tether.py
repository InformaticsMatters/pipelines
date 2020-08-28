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

import argparse, os, sys, json, traceback, time
import multiprocessing

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdForceFieldHelpers import *

from rdkit import RDLogger
from rdkit import rdBase

from pipelines.dimorphite.dimorphite_dl import run_with_mol_list


def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def enumerate_charges(mol, min_ph, max_ph):
    """
    Enumerate the charge states of this mol
    :param mol:
    :return:
    """

    mol_list = [mol]
    protonated_mols = run_with_mol_list(mol_list, min_ph=min_ph, max_ph=max_ph)

    return protonated_mols


chunk_size = 0
write_count = 0
file_count = 0
writer = None
embedding_failures_file = None


def get_writer(outfile_base):
    global writer
    global write_count
    global file_count

    if not writer:
        writer = Chem.SDWriter(outfile_base + '_' + f'{file_count:03}' + '.sdf')

    if write_count <= chunk_size:
        # print('  Using existing writer for', id, writer)
        return writer

    # we need to create a new writer
    # re-set the molecule count
    write_count = 0
    # close the old writer
    print('Closing writer', file_count)
    writer.close()
    file_count += 1

    writer = Chem.SDWriter(outfile_base + '_' + f'{file_count:03}' + '.sdf')
    # print('  Using new writer for', id, writer)
    return writer


def get_h_bonds(atom):
    l = []
    for b in atom.GetBonds():
        o = b.GetOtherAtom(atom)
        if o.GetAtomicNum() == 1:
            l.append(b.GetIdx())
    return tuple(l)


def get_num_h_neighbours(atom):
    count = 0
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != 1:
            count += 1
    return len(atom.GetNeighbors()) - count


def count_h_attachments_for_match(mol, match):
    res = []
    for m in match:
        atom = mol.GetAtomWithIdx(m)
        res.append(get_num_h_neighbours(atom))
    return tuple(res)


def count_h_attachements_for_mol(mol, match):
    molh = Chem.AddHs(mol)
    res = []
    for m in match:
        atom = molh.GetAtomWithIdx(m)
        if atom.GetAtomicNum() != 1:
            res.append(get_num_h_neighbours(atom))
    return tuple(res)


def count_attachments(mol, match):
    res = []
    for m in match:
        atom = mol.GetAtomWithIdx(m)
        res.append(get_h_bonds(atom))
    return tuple(res)


num_swap_success = 0
embedding_failures = 0

from random import randint


def do_embed(new_mol, coord_map, ret_vals):
    # for some strange reason you get the same conformer unless you specify a different seed.
    seed = randint(0, 100000000)
    ci = AllChem.EmbedMolecule(new_mol, coordMap=coord_map, ignoreSmoothingFailures=True, enforceChirality=False, randomSeed=seed)
    ret_vals['ci'] = ci
    ret_vals['mol'] = new_mol


def do_embed_with_timeout(new_mol, coord_map, timout_secs):
    manager = multiprocessing.Manager()
    ret_values = manager.dict()
    p = multiprocessing.Process(target=do_embed, name="EmbedMolecule", args=(new_mol, coord_map, ret_values))
    t0 = time.time()
    p.start()
    p.join(timout_secs)
    # If thread is active
    if p.is_alive():
        print("  Embedding is still running after", timout_secs, "seconds... let's kill it...")
        p.terminate()
        p.join()
        ci = -1
        mol = None
    else:
        t1 = time.time()
        ci = ret_values['ci']
        mol = ret_values['mol']
        mol.SetDoubleProp('EmbedTime', t1 - t0)

    return ci, mol


def multi_constrained_embed(mol, target, mcsQuery, getForceField=UFFGetMoleculeForceField,
                            num_conformers=1, timout_embed_secs=5, minimize_cycles=0):
    global num_swap_success
    global embedding_failures

    t0 = time.time()

    print(Chem.MolToSmiles(mol), Chem.MolToSmarts(mcsQuery))

    # In original constrained function, only a single match was found.
    # Here we look for multiple matches
    molMatches = mol.GetSubstructMatches(mcsQuery, uniquify=False)
    if not molMatches:
        raise ValueError("molecule doesn't match the MCS")
    targetMatches = target.GetSubstructMatches(mcsQuery, uniquify=True)
    if not targetMatches:
        raise ValueError("target doesn't match the MCS")
    print('  Found', len(molMatches), 'SSS matches to mol and', len(targetMatches), 'to target')

    # use the first target match for now
    targetMatch = targetMatches[0]
    print('  Target match:', targetMatch)

    mols_matches = []

    hs_for_core = count_h_attachements_for_mol(target, targetMatch)
    print('  COREHS', hs_for_core)

    # loop over our matches
    nhas = []

    index = 0
    for molMatch in molMatches:
        # print('  Match lengths:', len(targetMatch), len(molMatch))
        nha = count_h_attachments_for_match(mol, molMatch)
        atoms_to_switch = []
        if nha in nhas:
            print('  NHA', nha, 'SSS', molMatch, 'EXISTING')
            continue
        else:
            nhas.append(nha)
            for idx in range(len(hs_for_core)):
                hs_on_core_atoms = hs_for_core[idx]
                hs_on_matched_atoms = nha[idx]
                # check if the number of H atoms is different e.g. the atom is substituted compared to the core
                if hs_on_core_atoms != hs_on_matched_atoms:
                    # print('    ',idx, hs_on_core_atoms, hs_on_matched_atoms,c)
                    # Do we have at least one H on the atom - if so there is potential for swapping the substitutions
                    if hs_on_matched_atoms > 0:
                        atomidx = molMatch[idx]
                        atom = mol.GetAtomWithIdx(atomidx)
                        free_atoms = []
                        free_symbols = []
                        tethered_atoms = []
                        tethered_symbols = []
                        for other_atom in atom.GetNeighbors():
                            o_idx = other_atom.GetIdx()
                            o_sym = other_atom.GetSymbol() + str(other_atom.GetIdx())
                            if other_atom.GetIdx() in molMatch:
                                # print('  tethered', o_sym)
                                tethered_atoms.append(o_idx)
                                tethered_symbols.append(o_sym)
                            else:
                                # print('  free', o_sym)
                                free_atoms.append(o_idx)
                                free_symbols.append(o_sym)
                        print('  Atom:', atom.GetSymbol() + str(atomidx), 'Free:', free_symbols, 'Tethered:',
                              tethered_symbols)
                        # Include if there are 2 or more tethered atoms
                        # - otherwise it's free to rotate
                        if len(tethered_atoms) > 1:
                            # print('  GT1')
                            atoms_to_switch.append((idx, atomidx, molMatch[idx], hs_on_matched_atoms, tuple(free_atoms),
                                                    tuple(tethered_atoms)))

                        # Handle the special case of where there is only one tether but the bond is not rotatable
                        # e.g. a double bond, but not a methyl group where the 3 Hs are equivalent
                        if len(tethered_atoms) == 1:
                            # print('  EQ1')
                            bond = mol.GetBondBetweenAtoms(molMatch[idx], tethered_atoms[0])
                            if bond:  # should always be one, but just in case
                                # Don't know how to ask if a bond is rotatable.
                                # So instead ask if it is a double bond which won't be rotatable
                                if bond.GetBondType() == Chem.BondType.DOUBLE:
                                    print('  NON-ROTATABLE', bond.GetBondType(),
                                          atom.GetSymbol() + '->' + bond.GetOtherAtom(atom).GetSymbol())
                                    atoms_to_switch.append((idx, atomidx, molMatch[idx], hs_on_matched_atoms,
                                                            tuple(free_atoms), tuple(tethered_atoms)))

            print('  NHA', nha, 'SSS', molMatch, 'NEW', atoms_to_switch)

        free_counts = []
        for atom in mol.GetAtoms():
            free_count = 0
            for other_atom in atom.GetNeighbors():
                if other_atom.GetIdx() not in molMatch:
                    free_count += 1
            free_counts.append(free_count)

        mols_to_embed = enumerate_undefined_chirals(mol, free_counts)

        for new_mol in mols_to_embed:
            print('  Processing mol', Chem.MolToSmiles(new_mol))
            # print('  Target match:', symbolise_match(targetMatch, target))
            print('  Candid match:', symbolise_match(molMatch, new_mol))
            # print('  Candid atoms:', symbolise_mol(new_mol))
            coordMap = {}
            algMap = []
            targetConf = target.GetConformer(-1)
            for i, idxI in enumerate(targetMatch):
                # coordMap[idxI] = coreConf.GetAtomPosition(i)
                coordMap[molMatch[i]] = targetConf.GetAtomPosition(idxI)
                algMap.append((molMatch[i], targetMatch[i]))
            # print('  CoordMap:', dump_coord_map(coordMap))

            try:
                conf_count = 0
                for conf in range(num_conformers):
                    clone = Chem.RWMol(new_mol)
                    ci, embed_mol = do_embed_with_timeout(clone, coordMap, timout_embed_secs)

                    if ci < 0:
                        print('  WARNING: Could not embed molecule.')
                        embedding_failures += 1
                    else:
                        # rotate the embedded conformation onto the core:
                        rms = AlignMol(embed_mol, target, atomMap=algMap)

                        # process the original embedded molecule
                        if minimize_cycles > 0:
                            minimize_mol(embed_mol, target, molMatch, targetMatch, algMap, getForceField, minimize_cycles)
                        new_mol.SetIntProp('Index', index)
                        new_mol.SetIntProp('Conformer', conf)
                        mols_matches.append((embed_mol, molMatch))
                        conf_count += 1
                print('   ', conf_count, 'conformers generated')
            except:
                embedding_failures += 1
                print('  ERROR: Failed to Embed molecule')
                traceback.print_exc()

        if not mols_matches:
            # all embedding failed. Maybe due to troublesome stereo with highly constrained structures
            # so let's do one last attempt without the stereo
            ci, embed_mol = do_embed_with_timeout(Chem.RWMol(mol), coordMap, timout_embed_secs)
            if ci < 0:
                print('WARNING: Could not embed molecule.')
                embedding_failures += 1
            else:
                # rotate the embedded conformation onto the core:
                rms = AlignMol(embed_mol, target, atomMap=algMap)

                # process the original embedded molecule
                if minimize_cycles > 0:
                    minimize_mol(embed_mol, target, molMatch, targetMatch, algMap, getForceField, minimize_cycles)
                embed_mol.SetIntProp('Index', index)
                embed_mol.SetIntProp('Conformer', 0)
                mols_matches.append((embed_mol, molMatch))

        index += 1

    t1 = time.time()
    # print('  Embedding took:', t1 - t0)
    # Return a list of tuples of (mol, match)
    return mols_matches


def dump_coord_map(coordMap):
    res = ''
    for idx, point in coordMap.items():
        res += str(idx) + '->[' + str(point.x) + ',' + str(point.y) + ',' + str(point.z) + '] '
    return res

def symbolise_match(match, mol):
    r = []
    for i in match:
        r.append(mol.GetAtomWithIdx(i).GetSymbol() + str(i))
    return r

def symbolise_mol(mol):
    r = []
    for i in range(mol.GetNumAtoms()):
        r.append(mol.GetAtomWithIdx(i).GetSymbol() + str(i))
    return r

def enumerate_undefined_chirals(mol, free_counts):
    """
    Enumerate undefined chiral centres that are points of substitution.
    Chiral centres that are already specifically defined are left untouched.
    :param mol: The mol to enumerate
    :param atoms_to_switch: The points of substitution
    :return:
    """

    global num_swap_success

    chiral_targets = []
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    print('  Chiral centers:', chiral_centers, 'Free atom counts:', free_counts)

    for i, cc in enumerate(chiral_centers):
        if cc[1] == '?':  # only enumerate if the stereo is undefined
            if free_counts[i] > 1:  # only enumerate if there are at least 2 free atoms
                chiral_targets.append(cc[0])

    if chiral_targets:
        mols_to_embed = []
        mol_cw = Chem.RWMol(mol)
        mol_cw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        mols_to_embed.append(mol_cw)
        mol_ccw = Chem.RWMol(mol)
        mol_ccw.GetAtomWithIdx(chiral_targets[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        mols_to_embed.append(mol_ccw)
        for j in range(1, len(chiral_targets)):
            new_enum_mols = []
            for m in mols_to_embed:
                nm1 = Chem.RWMol(m)
                nm1.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                new_enum_mols.append(nm1)
                nm2 = Chem.RWMol(m)
                nm2.GetAtomWithIdx(chiral_targets[j]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                new_enum_mols.append(nm2)
                mols_to_embed = new_enum_mols
                num_swap_success += len(mols_to_embed)
    else:
        mols_to_embed = [Chem.RWMol(mol)]
    print('  Mols to embed:', len(mols_to_embed))
    return mols_to_embed


def minimize_mol(mol, target, molMatch, targetMatch, algMap, getForceField, cycles):
    ff = getForceField(mol, confId=-1)
    conf = target.GetConformer()
    for i, atIdx in enumerate(targetMatch):
        p = conf.GetAtomPosition(atIdx)
        pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
        ff.AddDistanceConstraint(pIdx, molMatch[i], 0, 0, 100.0)
    ff.Initialize()
    n = cycles
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        n -= 1
    # realign
    rms = AlignMol(mol, target, atomMap=algMap)
    mol.SetDoubleProp('EmbedRMS', rms)


def find_best_mcs(hit, mol):
    best_score = 0
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=True, ringMatchesRingOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 1
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=True, ringMatchesRingOnly=False,
                         atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 2
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=False, ringMatchesRingOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 3
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=False, ringMatchesRingOnly=False,
                         atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrder)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 4

    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=True, ringMatchesRingOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 5
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=True, ringMatchesRingOnly=False,
                         atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 6
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=False, ringMatchesRingOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 7
    mcs = rdFMCS.FindMCS([hit, mol], completeRingsOnly=False, ringMatchesRingOnly=False,
                         atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny)
    score = score_mcs(hit, mol, mcs)
    if score > best_score:
        best_score = score
        best_mcs = mcs
        best_index = 8
    print("  Best MCS is", best_index, 'with score', best_score)
    return best_mcs


def score_mcs(hit, mol, mcs):
    smartsMol = Chem.MolFromSmarts(mcs.smartsString)
    hitMatches = hit.GetSubstructMatch(smartsMol)
    molMatches = mol.GetSubstructMatch(smartsMol)
    # print('  Matches:', hitMatches, molMatches)
    score = 0
    for hitMatch, molMatch in zip(hitMatches, molMatches):
        score += 1
        if hit.GetAtomWithIdx(hitMatch).GetAtomicNum() == mol.GetAtomWithIdx(molMatch).GetAtomicNum():
            score += 1
    if mol.GetAtomWithIdx(molMatch).IsInRing():
        score += 1
    print('  Score:', score)
    return score


def execute(smi, hit_molfile, outfile_base, min_ph=None, max_ph=None, max_inputs=0, max_outputs=0, modulus=0,
            num_conformers=1, timout_embed_secs=5,
            complete_rings_only=True, ring_matches_ring_only=True,
            atom_compare=rdFMCS.AtomCompare.CompareElements, bond_compare=rdFMCS.BondCompare.CompareOrder,
            minimize_cycles=0
            ):
    global write_count

    GetFF = lambda x, confId=-1: AllChem.MMFFGetMoleculeForceField(x, AllChem.MMFFGetMoleculeProperties(x),
                                                                   confId=confId)

    hit = Chem.MolFromMolFile(hit_molfile)

    num_mols = 0
    num_processed = 0
    num_outputs = 0
    num_errors = 0
    complete_failures = 0

    with open(smi, 'r') as fsmi:

        for line in fsmi:
            t0 = time.time()

            smiles = line.strip()

            if 0 < modulus and num_mols % modulus != 0:
                num_mols += 1
                continue

            if 0 < max_inputs <= num_mols:
                break

            if 0 < max_outputs <= num_processed:
                break

            num_mols += 1
            num_processed += 1
            num_added = 0

            w = get_writer(outfile_base)

            print('Processing', num_mols, num_processed, smiles, Chem.MolToSmiles(hit))

            try:
                mol = Chem.MolFromSmiles(smiles)

                if min_ph and max_ph:
                    enumerated_mols = enumerate_charges(mol, min_ph, max_ph)
                    print('  Enumerated', len(enumerated_mols), [Chem.MolToSmiles(x) for x in enumerated_mols])
                else:
                    enumerated_mols = [mol]

                mcs0 = rdFMCS.FindMCS([hit, mol], completeRingsOnly=complete_rings_only,
                                      ringMatchesRingOnly=ring_matches_ring_only,
                                      atomCompare=atom_compare, bondCompare=bond_compare)
                # mcs0 = rdFMCS.FindMCS([hit, mol], completeRingsOnly=False, ringMatchesRingOnly=False)
                # score_mcs(hit, mol, mcs0)

                # mcs0 = find_best_mcs(hit, mol)

                mcsQuery = Chem.MolFromSmarts(mcs0.smartsString)

                count = 0
                for ligand in enumerated_mols:
                    molh = Chem.AddHs(ligand)
                    mol_match_tuple = multi_constrained_embed(molh, hit, mcsQuery, num_conformers=num_conformers,
                                                              timout_embed_secs=timout_embed_secs, minimize_cycles=minimize_cycles)
                    print(' ', len(mol_match_tuple), 'mols tethered')

                    for t_mol, match in mol_match_tuple:
                        count += 1
                        t_mol.SetProp('_Name', smiles)
                        # t_mol.SetIntProp('Index', count)
                        if hit.HasProp('_Name'):
                            t_mol.SetProp('Hit', hit.GetProp('_Name'))
                        else:
                            t_mol.SetProp('Hit', hit_molfile[:-4])
                        tethers = ','.join([str(x + 1) for x in match])
                        t_mol.SetProp('TETHERED ATOMS', tethers)
                        # print(' Tethers: ', tethers)

                        w.write(t_mol)
                        write_count += 1
                        num_outputs += 1
                        num_added += 1

            except Exception as e:
                num_errors += 1
                print('  Error: ', smiles, e)
                traceback.print_exc()

            t1 = time.time()
            print('  processed mol', num_processed, 'in', t1 - t0, 'secs. Generated', num_added, 'mols')
            if num_added == 0:
                embedding_failures_file.write(line)
                embedding_failures_file.flush()
                complete_failures += 1

        w.close()

        print('Totals - inputs:', num_mols, 'processed:', num_processed, 'outputs:', num_outputs, 'chiral swaps:',
              num_swap_success, 'errors:', num_errors, 'embed failures:', embedding_failures, 'failures:',
              complete_failures, 'files:', file_count + 1)


# atom_compares = {
#     'elements': rdFMCS.AtomCompare.CompareElements,
#     'any': rdFMCS.AtomCompare.CompareAny,
#     'anyHeavyAtom': rdFMCS.AtomCompare.CompareAnyHeavyAtom,
#     'isotopes': rdFMCS.AtomCompare.CompareIsotopes
# }
#
# bond_copares = {
#     'any': rdFMCS.BondCompare.CompareAny,
#     'order': rdFMCS.BondCompare.CompareOrder,
#     'exact': rdFMCS.BondCompare.CompareExact
# }

def main():
    """
    Example usage:
    python -m pipelines.xchem.prepare_tether --smi ../../data/mpro/Mpro-x0387_0.smi --mol ../../data/mpro/Mpro-x0387_0.mol -o TETHERED --max-inputs 500 --chunk-size 100

    :return:
    """

    global chunk_size
    global embedding_failures_file

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)
    print('RDKit version:', rdBase.rdkitVersion)

    parser = argparse.ArgumentParser(description='Tether prep - prepare candidates for docking')

    parser.add_argument('--smi', help='SMILES containing the expanded candidates for a hit)')
    parser.add_argument('--mol', help='Molfile containing the hit to tether to)')
    parser.add_argument('-o', '--outfile', default='Tethered',
                        help='Base name for results SDF file (will generate something like Tethered_Mpro-x0072_000.sdf)')
    parser.add_argument('--min-ph', type=float, help='The min pH to consider')
    parser.add_argument('--max-ph', type=float, help='The max pH to consider')
    parser.add_argument('-c', '--chunk-size', type=int, default=200, help='Chunk size for files')
    parser.add_argument('--max-inputs', type=int, default=0, help='Max number of molecules to process')
    parser.add_argument('--max-outputs', type=int, default=0, help='Max number of records to output')
    parser.add_argument('--modulus', type=int, default=0, help='Process only mols with this modulus')
    parser.add_argument('--num-conformers', type=int, default=1,
                        help='Generate this number of conformers for each tethering')
    parser.add_argument('--timeout-embed', type=int, default=5, help='Timeout in seconds to apply to limit embedding')
    parser.add_argument('--ring-matches-ring-only', action="store_true", help='Set ringMatchesRingOnly MCS property to True')
    parser.add_argument('--complete-rings-only', action="store_true", help='Set completeRingsOnly MCS property to True')
    parser.add_argument('--atom-compare', default='CompareElements', help='atomCompare MCS property',
                        choices=['CompareAny', 'CompareAnyHeavyAtom', 'CompareElements', 'CompareIsotopes'])
    parser.add_argument('--bond-compare', default='CompareOrder', help='bondCompare MCS property',
                        choices=['CompareAny', 'CompareOrder', 'CompareOrderExact'])
    parser.add_argument("--minimize", type=int, default=0, help="number of minimisation cycles")


    args = parser.parse_args()
    log("Tether prep args: ", args)

    chunk_size = args.chunk_size

    min_ph = args.min_ph
    max_ph = args.max_ph
    smi = args.smi
    mol = args.mol
    outfile = args.outfile
    max_inputs = args.max_inputs
    max_outputs = args.max_outputs
    modulus = args.modulus
    num_conformers = args.num_conformers
    timout_embed_secs = args.timeout_embed
    ring_matches_ring_only = args.ring_matches_ring_only
    complete_rings_only = args.complete_rings_only
    atom_compare = rdFMCS.AtomCompare.names[args.atom_compare]
    bond_compare = rdFMCS.BondCompare.names[args.bond_compare]
    minimize_cycles = args.minimize

    embedding_failures_file = open(outfile + '_embedding_failures.smi', 'w')

    # Dimporphite needs to use argparse with its own arguments, not messed up with our arguments
    # so we store the original args
    orig_sys_argv = sys.argv[:]

    # Remove all the parameters, keeping only the filename (first one) so that
    # dimorphite is unaware of any previous commandline parameters.
    sys.argv = sys.argv[:1]

    execute(smi, mol, outfile, min_ph=min_ph, max_ph=max_ph,
            max_inputs=max_inputs, max_outputs=max_outputs, modulus=modulus,
            num_conformers=num_conformers, timout_embed_secs=timout_embed_secs,
            complete_rings_only=complete_rings_only, ring_matches_ring_only=ring_matches_ring_only,
            atom_compare=atom_compare, bond_compare=bond_compare,
            minimize_cycles=minimize_cycles)

    embedding_failures_file.close()

    print('Finished')


if __name__ == "__main__":
    main()
