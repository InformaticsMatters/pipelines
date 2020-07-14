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

import os, argparse

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from datetime import datetime

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

field_XCosRefMols = "XCos_RefMols"
field_XCosNumHits = "XCos_NumHits"
field_XCosScore1 = "XCos_Score1"
field_XCosScore2 = "XCos_Score2"
field_XCosScore3 = "XCos_Score3"


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

def getNumberfeats(mol):

    featLists = []
    rawFeats = fdef.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're intereted in
    featLists.append([f for f in rawFeats if f.GetFamily() in keep])

    return len(featLists)


def getFeatureMapXCOS(mol_list):
    allFeats = []
    for m in mol_list:

        rawFeats = fdef.GetFeaturesForMol(m)
        featDeats = [(f.GetType(),
                      f.GetPos().x,
                      f.GetPos().y,
                      f.GetPos().z) for f in rawFeats if f.GetFamily() in keep]

        allFeats.append(featDeats)


    feature_map_df = pd.DataFrame([t for lst in allFeats for t in lst],
                                  columns =['featType', 'x', 'y', 'z'])

    return feature_map_df


def getFeatureAgg(feature_map_df, rad_thresh):

    # Group data into unique feature types
    grouped_df = feature_map_df.groupby('featType')

    data_to_add = []

    for group_name, df_group in grouped_df:

        # Reset index df
        df_group = df_group.reset_index()

        if len(df_group) == 1:

            data_to_add.append(df_group)

        if len(df_group) > 1:

            # Get feature name
            feat_name = df_group.featType.unique()[0]

            # Use radius neighbours to find features within
            # spere with radius thresh
            neigh = NearestNeighbors(radius=rad_thresh)

            while len(df_group) > 0:

                neigh.fit(df_group[['x','y','z']])

                # Get distances and indices of neigbours within radius threshold
                rng = neigh.radius_neighbors()
                neigh_dist = rng[0][0]
                neigh_indices = rng[1][0]

                # Append the first index - NB clustering done relative to index 0
                neigh_indices = list(np.append(0, neigh_indices))

                # Calculate average x,y,z coords for features in similar loc
                x_avg = np.mean(df_group.iloc[neigh_indices].x)
                y_avg = np.mean(df_group.iloc[neigh_indices].y)
                z_avg = np.mean(df_group.iloc[neigh_indices].z)

                # Add feature with average x, y and z values
                new_row = [(feat_name, x_avg, y_avg, z_avg)]

                cluster_df = pd.DataFrame(data=new_row, columns = ['featType', 'x', 'y', 'z'])

                data_to_add.append(cluster_df)

                # Remove indices of clustered neigbours
                df_group = df_group.drop(df_group.index[neigh_indices])

    # Create single DF from list of dfs
    clustered_df = pd.concat(data_to_add)

    return clustered_df


# This is the main XCOS function
def getReverseScores(clustered_df, mols, frags, no_clustered_feats, rad_threshold, COS_threshold, writer):

    for mol in mols:

        # Get the bits
        compound_bits = getBits(mol)

        # We are going to include a feature mapping score, where the
        # number of features of the compound matching the clustered feats
        # within a threshold are found

        # Get feature map of compound bits as df
        feature_map_bits = getFeatureMapXCOS(compound_bits)

        # Group data into unique feature types
        grouped_df = feature_map_bits.groupby('featType')

        no_feats_matched = []
        dist_feats_matched = []

        # Use radius neighbours to find features within
        # sphere with radius thresh
        neigh = NearestNeighbors(radius=rad_threshold)

        # Loop through grouped features
        for group_name, df_group in grouped_df:

            # Get feat name
            feat_name = df_group.featType.unique()[0]

            # Get similar feats from cluster df
            cluster_test = clustered_df[clustered_df.featType == feat_name]

            # Reset index df
            df_group = df_group.reset_index()

            if len(cluster_test) == 1:

                # Calculate distances
                x1_sub_x2 = (cluster_test.iloc[0].x - df_group.iloc[0].x) ** 2
                y1_sub_y2 = (cluster_test.iloc[0].y - df_group.iloc[0].y) ** 2
                z1_sub_z2 = (cluster_test.iloc[0].z - df_group.iloc[0].z) ** 2

                diff_sum = x1_sub_x2 + y1_sub_y2 + z1_sub_z2

                dist = diff_sum ** 0.5

                if dist < rad_threshold:
                    # Let's get the number of feats matched
                    no_feats_matched.append(1)

                    # Let's get the distance of the feats matched
                    dist_feats_matched.append([dist])

            if len(cluster_test) > 1:
                neigh.fit(cluster_test[['x', 'y', 'z']])

                while len(df_group) > 0:
                    # Get distances and indices of neigbours within radius threshold
                    feat_coords = [[df_group.iloc[0].x, df_group.iloc[0].y, df_group.iloc[0].z]]
                    rng = neigh.radius_neighbors(feat_coords)

                    neigh_dist = rng[0][0]
                    neigh_indices = rng[1][0]

                    # Let's get the number of feats matched
                    no_feats_matched.append(len(neigh_indices))

                    # Remove index 0 of df_group
                    df_group = df_group.drop(df_group.index[0])

        # Get total number of feat matches
        no_feats = np.sum(no_feats_matched)

        all_scores = []

        for bit in compound_bits:

            # Get number of bit atoms
            no_bit_atoms = bit.GetNumAtoms()

            scores = []
            for frag_mol in frags:
                # NB reverse SuCOS scoring
                fm_score = getFeatureMapScore(bit, frag_mol)
                fm_score = np.clip(fm_score, 0, 1)
                protrude_dist = rdShapeHelpers.ShapeProtrudeDist(bit, frag_mol, allowReordering=False)
                protrude_dist = np.clip(protrude_dist, 0, 1)

                reverse_SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

                # Get number of feats from bit for scaling score
                no_bit_feats = getNumberfeats(bit)

                # Get some info and append to list
                frag_name = frag_mol.GetProp('_Name').strip('Mpro-')

                scores.append((frag_name, reverse_SuCOS_score, no_bit_atoms, no_bit_feats))

            all_scores.append(scores)

            list_dfs = []
            for score in all_scores:
                df = pd.DataFrame(data=score, columns=['Fragment', 'Score', 'No_bit_atoms', 'No_bit_feats'])
                # Get maximum scoring fragment for bit match
                df = df[df['Score'] == df['Score'].max()]
                list_dfs.append(df)

            final_df = pd.concat(list_dfs)

            # Get total bit score and some denominator terms
            bits_score = (final_df.No_bit_atoms * final_df.Score).sum()
            total_atoms = final_df.No_bit_atoms.sum()
            feat_match_fraction = no_feats / no_clustered_feats

            # Score 1: the score is scaled by the number of bit atoms
            score_1 = bits_score

            # Score 2: the score is scaled by the number of bit atoms
            # penalised by the fraction of feats matched
            # the to total number feats clustered
            score_2 = score_1 * feat_match_fraction

            # Score 3: the score is determined by the fraction of matching
            # features to the clustered features within a threshold. This
            # should yield similar values to Tim's Featurestein method?
            score_3 = feat_match_fraction

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
        mol.SetProp(field_XCosScore2, "{:.4f}".format(score_2))
        mol.SetProp(field_XCosScore3, "{:.4f}".format(score_3))

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

    feature_map_df =  getFeatureMapXCOS(frag_mol_list)
    utils.log('Feature map dataframe shape:', feature_map_df.shape)

    # Set radius threshold
    rad_thresh = 1.5

    # Aggregate features using nearest neigbours algo
    clustered_df = getFeatureAgg(feature_map_df, rad_thresh=rad_thresh)
    utils.log('Clustered dataframe shape:', clustered_df.shape)
    no_clustered_feats = len(clustered_df)

    #clustered_df, mols, rad_threshold, COS_threshold, writer
    getReverseScores(clustered_df, molecules, frag_mol_list, no_clustered_feats, 1.0, 0.50, writer)


def main():

    # Example usage:
    #  python -m pipelines.xchem.xcos -f ../../data/mpro/hits-17.sdf.gz -i ../../data/mpro/poses.sdf.gz  -o xcos

    global fmaps

    parser = argparse.ArgumentParser(description='XCos scoring with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-f', '--fragments', required=True, help='Fragments to compare')
    parser.add_argument('--fragments-format', help='Fragments format')
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
    clsMappings[field_XCosScore2] = "java.lang.Float"
    clsMappings[field_XCosScore3] = "java.lang.Float"
    fieldMetaProps.append({"fieldName":field_XCosRefMols,   "values": {"source":source, "description":"XCos reference fragments"}})
    fieldMetaProps.append({"fieldName":field_XCosNumHits,   "values": {"source":source, "description":"XCos number of hits"}})
    fieldMetaProps.append({"fieldName":field_XCosScore1,   "values": {"source":source, "description":"XCos score 1"}})
    fieldMetaProps.append({"fieldName":field_XCosScore2,   "values": {"source":source, "description":"XCos score 2"}})
    fieldMetaProps.append({"fieldName":field_XCosScore3,   "values": {"source":source, "description":"XCos score 3"}})


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
