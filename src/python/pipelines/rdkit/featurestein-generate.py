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
Ligand pose scoring using 'FeatureStein'.
This module generates a merged feature map from a set of 3D ligands.
The output is a pickle of the merged feature map that can be read by the featurestein-score.py module to
generate scores.
"""

from __future__ import print_function
import argparse, os, sys, gzip, pickle

from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.FeatMaps.FeatMapUtils import CombineFeatMaps

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


### start function definitions #########################################

ffact = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

fmParams = {}
for k in ffact.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

exclude = ()

def filterFeature(f):
    if f.GetFamily() in exclude:
        return None
    else:
        return f

def getRawFeatures(mol):
    rawFeats = ffact.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filterFeature, rawFeats))
    return filtered

def getFeatureMap(mol):
    feats = getRawFeatures(mol)
    return FeatMaps.FeatMap(feats=feats, weights=[1]*len(feats),params=fmParams)

def score_featmaps(fm1, fm2):
    "Generate the score for 2 feature maps"
    return fm1.ScoreFeats(fm2.GetFeatures()) / fm1.GetNumFeatures()

def build_feat_data(mols):
    "Build the feature maps and do the all vs. all comparison"
    fmaps = []
    scores = []
    for mol1 in mols:
        fm1 = getFeatureMap(mol1)
        fmaps.append(fm1)
        row = []
        for mol2 in mols:
            fm2 = getFeatureMap(mol2)
            score = score_featmaps(fm1, fm2)
            row.append(score)
            #print(len(data), len(row), score)
        scores.append(row)
    return fmaps, scores

def find_closest(scores):
    #print('Find closest for', len(scores), len(scores[0]))
    best_score = 0
    for i in range(len(scores)):
        for j in range(len(scores)):
            if i == j:
                continue
            score = scores[i][j]
            if score > best_score:
                best_score = score
                best_row = i
                best_col = j
    return best_score, best_row, best_col

def merge_feat_maps(fmaps, scores):
    "Merge the 2 closest feature maps, remove them form the data and replace with the merged feature map"
    best_score, best_row, best_col = find_closest(scores)
    #print(best_score, best_row, best_col)
    feat1 = fmaps[best_row]
    feat2 = fmaps[best_col]
    utils.log('Merging', best_row, 'and', best_col, 'with score', best_score, '#features:', feat1.GetNumFeatures(), feat2.GetNumFeatures())
    merged = CombineFeatMaps(feat1, feat2, mergeMetric=1, mergeTol=1.5, dirMergeMode=0)
    # need to make sure we delete the biggest index first to avoid changing the smaller index
    if best_row > best_col:
        a = best_row
        b = best_col
    else:
        a = best_col
        b = best_row

    #print('Initial:', len(fmaps), len(scores), ','.join([str(len(x)) for x in scores]))
    del fmaps[a]
    del fmaps[b]
    del scores[a]
    del scores[b]
    for row in scores:
        del row[a]
        del row[b]

    merged_scores = []
    for i in range(len(fmaps)):
        fmap = fmaps[i]
        score1 = score_featmaps(fmap, merged)
        score2 = score_featmaps(merged, fmap)
        scores[i].append(score1)
        merged_scores.append(score2)

    fmaps.append(merged)
    merged_scores.append(score_featmaps(merged, merged))
    scores.append(merged_scores)


def process(inputs, fname):

    mols = [m for m in inputs if m]
    fmaps, scores = build_feat_data(mols)
    merged_fmaps = fmaps.copy()
    utils.log('Processing', len(fmaps), 'molecules')
    while len(merged_fmaps) > 1:
        merge_feat_maps(merged_fmaps, scores)
    merged_fmap = merged_fmaps[0]
    pickle.dump(merged_fmap, open(fname, "wb" ))
    utils.log('Wrote merged feature map with', merged_fmap.GetNumFeatures(), 'features as pickle to', fname)

    return len(mols), merged_fmap.GetNumFeatures()

### start main execution #########################################

def main():

    global fmaps

    parser = argparse.ArgumentParser(description='FeatureStein generation with RDKit')
    parameter_utils.add_default_input_args(parser)
    parser.add_argument('-f', '--feat-map', default='featurestein.p', help='Name of pickle to generate')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')

    args = parser.parse_args()
    utils.log("FeatureStein Args: ", args)

    inputs_file, inputs_supplr = rdkit_utils. \
        default_open_input(args.input, args.informat)

    # this does the processing
    num_mols, num_feats = process(inputs_supplr, args.feat_map)

    inputs_file.close()

    if args.metrics:
        utils.write_metrics(output_base, {'__StatusMessage__': 'Generated ' + num_feats + ' from ' + num_mols + ' molecules',
                                          '__InputCount__':num_mols, 'RDKitFeatureMap':num_mols})


if __name__ == "__main__":
    main()
