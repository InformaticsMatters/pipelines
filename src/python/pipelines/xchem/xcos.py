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
XCos fragments a molecule and compares the fragments to a series of ligands.
The XCos code was written by Warren Thompson <warren.thompson@diamond.ac.uk>
"""

from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem import BRICS
from rdkit import Chem, rdBase
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit import RDConfig
from rdkit.Chem import rdFMCS

import os, argparse

import numpy as np
import pandas as pd

from datetime import datetime

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

field_XCosRefMols = "XCos_RefMols"
field_XCosNumHits = "XCos_NumHits"
field_XCosScore1 = "XCos_Score1"


date = datetime.today().strftime('%Y-%m-%d')

def getBits(mol):
    '''

    Parameters
    ----------
    mol : rdkit mol object to be broken up into fragments by breaking
    rotable bonds

    Returns
    -------
    mols : A list of rdkit mol objects

    '''
    # find the rotatable bonds
    bonds = mol.GetSubstructMatches(RotatableBondSmarts)

    bonds = [((x,y),(0,0)) for x,y in bonds]
    p = BRICS.BreakBRICSBonds(mol,bonds=bonds)

    mols = [mol for mol in Chem.GetMolFrags(p,asMols=True)]

    return mols

# We need to start by building a FeatureFactory object which defines
# the set of pharmacophore features being used.
# We'll use this to find features on the molecules.
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,
                                                'BaseFeatures.fdef'))


# Set default paramters for selecting points in feature map
fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

# List of feature families that we want to use
keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')


def getFeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):
    try:
        featLists = []
        for m in [small_m, large_m]:
            rawFeats = fdef.GetFeaturesForMol(m)
            # filter that list down to only include the ones we're intereted in
            featLists.append([f for f in rawFeats if f.GetFamily() in keep])
        fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
        fms[0].scoreMode = score_mode
        fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
        return fm_score
    except ZeroDivisionError:
        return 0


# This is the main XCOS function
def getReverseScores(mols, frags, COS_threshold, writer):

    for mol in mols:

        # Get the bits
        compound_bits = getBits(mol)

        all_scores = []

        for bit in compound_bits:

            # Get number of bit atoms
            no_bit_atoms = bit.GetNumAtoms()

            scores = []

            for frag_mol in frags:

                # Get frag name for linking to score
                frag_name = frag_mol.GetProp('_Name').strip('Mpro-')
                
                # Check if MCS yield > 0 atoms
                mcs_match = rdFMCS.FindMCS([bit,frag_mol],ringMatchesRingOnly=True,matchValences=True)
                
                # Get number of atoms in MCS match found
                no_mcs_atoms = Chem.MolFromSmarts(mcs_match.smartsString).GetNumAtoms()

                if no_mcs_atoms == 0:

                    scores.append((frag_name, 0, no_bit_atoms))
                
                if no_mcs_atoms > 0:

                    # NB reverse SuCOS scoring
                    fm_score = getFeatureMapScore(bit, frag_mol)
                    fm_score = np.clip(fm_score, 0, 1)

                    # Change van der Waals radius scale for stricter overlay                     
                    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(bit, frag_mol,
                                                                     allowReordering=False,
                                                                     vdwScale=0.2)
                    protrude_dist = np.clip(protrude_dist, 0, 1)

                    reverse_SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

                    scores.append((frag_name, reverse_SuCOS_score, no_bit_atoms))

            all_scores.append(scores)

            list_dfs = []

            for score in all_scores:

                df = pd.DataFrame(data=score, columns=['Fragment', 'Score', 'No_bit_atoms'])
                
                # Get maximum scoring fragment for bit match
                df = df[df['Score'] == df['Score'].max()]
                list_dfs.append(df)

            final_df = pd.concat(list_dfs)

            # Score 1: the score is scaled by the number of bit atoms
            score_1 = (final_df.No_bit_atoms * final_df.Score).sum()
           
            # Let's only get frags above a threshold
            final_df = final_df[final_df.Score > COS_threshold]

            # Let#s sort the df by increasing score
            final_df = final_df.sort_values(by=['Score'], ascending=False)

            # Get the unique fragments above threshold
            all_frags = pd.unique(final_df.Fragment)

        # Add props we want
        mol.SetProp(field_XCosRefMols, ','.join(all_frags))
        mol.SetIntProp(field_XCosNumHits, len(all_frags))
        mol.SetProp(field_XCosScore1, "{:.4f}".format(score_1))

        # Write to file
        writer.write(mol)
        writer.flush()

def process(molecules, fragments, writer):

    frag_mol_list = []
    errors = 0
    for m in fragments:
        if m is None:
            errors += 1
        else:
            frag_mol_list.append(m)
    if errors:
        utils.log(errors, 'molecules failed to load. Using', len(frag_mol_list), 'fragments.')
    else:
        utils.log('Using', len(frag_mol_list), 'fragments. No errors')

    #mols, frags, COS_threshold, writer
    getReverseScores(molecules, frag_mol_list, 0.40, writer)


def main():

    # Example usage:
    #  python -m pipelines.xchem.xcos -f ../../data/mpro/hits-17.sdf.gz -i ../../data/mpro/poses.sdf.gz  -o xcos

    global fmaps

    parser = argparse.ArgumentParser(description='XCos scoring with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-f', '--fragments', required=True, help='Fragments to compare')
    parser.add_argument('-ff', '--fragments-format', help='Fragments format')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')

    args = parser.parse_args()
    utils.log("XCos Args: ", args)

    source = "xcos.py"
    datasetMetaProps = {"source":source, "description": "XCos scoring using RDKit " + rdBase.rdkitVersion}

    clsMappings = {}
    fieldMetaProps = []

    clsMappings[field_XCosRefMols] = "java.lang.String"
    clsMappings[field_XCosNumHits] = "java.lang.Integer"
    clsMappings[field_XCosScore1] = "java.lang.Float"

    fieldMetaProps.append({"fieldName":field_XCosRefMols,   "values": {"source":source, "description":"XCos reference fragments"}})
    fieldMetaProps.append({"fieldName":field_XCosNumHits,   "values": {"source":source, "description":"XCos number of hits"}})
    fieldMetaProps.append({"fieldName":field_XCosScore1,   "values": {"source":source, "description":"XCos score 1"}})
    
    frags_input,frags_suppl = rdkit_utils.default_open_input(args.fragments, args.fragments_format)

    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                                                                  'xcos', args.outformat,
                                                                  valueClassMappings=clsMappings,
                                                                  datasetMetaProps=datasetMetaProps,
                                                                  fieldMetaProps=fieldMetaProps,
                                                                  compress=not args.no_gzip)

    # this does the processing
    process(inputs_supplr, frags_suppl, writer)

    writer.close()


if __name__ == "__main__":
    main()
