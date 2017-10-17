#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
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

import logging
from rdkit import Chem
from rdkit.Chem import Descriptors
from pipelines.utils import utils

def fragment(mol, mode, quiet=False):
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol
    else:
        # TODO - handle ties
        biggest_index = -1
        i = 0
        if mode == 'hac':
            biggest_count = 0
            for frag in frags:
                hac = frag.GetNumHeavyAtoms()
                if hac > biggest_count:
                    biggest_count = hac
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            if not quiet:
                utils.log("Chose fragment", biggest_index, "from", len(frags), "based on HAC")
        elif mode == 'mw':
            biggest_mw = 0
            for frag in frags:
                mw = Descriptors.MolWt(frag)
                if mw > biggest_mw:
                    biggest_mw = mw
                    biggest_mol = frag
                    biggest_index = i
                i+=1
            if not quiet:
                utils.log("Chose fragment", biggest_index, "from", len(frags), "based on MW")
        else:
            raise ValueError('Invalid fragment mode:',mode)

        # copy the properties across
        for name in mol.GetPropNames():
            biggest_mol.SetProp(name, mol.GetProp(name))
        name = mol.GetProp("_Name")
        if name:
            biggest_mol.SetProp("_Name", name)
        return biggest_mol

def fragmentAndFingerprint(suppl, mols, fps, descriptor, fragmentMethod='hac', outputFragment=False, quiet=False):
    """
    Fragment the molecule if it has multiple fragments and generate fingerprints on the fragment.

    :param suppl: MolSupplier from which to read the molecules
    :param mols: List to which the molecules are added
    :param fps: List to which the fingerprints are added
    :param descriptor: Function to generate the fingerprints from the molecule
    :param fragmentMethod: The fragmentation method to use when there are multiple fragments (hac or mw)
    :param outputFragment: Boolean that specifies whether to write the fragment or the original molecule to the mols list
    :param quiet: Quiet mode
    :return: The number of errors encountered
    """
    errors = 0
    for mol in suppl:
        try:
            if mol is not None:
                frag = fragment(mol, fragmentMethod, quiet)
                d = descriptor(frag)
                if d:
                    if outputFragment:
                        mols.append(frag)
                    else:
                        mols.append(mol)
                    fps.append(d)
                    continue
        except:
            logging.exception('')
        errors += 1
    return errors