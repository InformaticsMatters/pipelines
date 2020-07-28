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
FeatureStein is a merged RDKit feature map that estimates the overlap of a ligand with a set of ligands (e.g. fragment
screening hits) based on the RDKit feature maps.
This modules generates merged features maps for a set of fragments and then scores a a set of ligands with that merged feature map.
See featurestein_generate.py and featurestein_generate.py for how to separate these into separate processes.
"""

from __future__ import print_function
import argparse, os, sys, gzip, pickle, traceback

from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.FeatMaps.FeatMapUtils import CombineFeatMaps

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


### start function definitions #########################################

field_FeatureSteinQualityScore = "FeatureStein_Qual"
field_FeatureSteinQuantityScore = "FeatureStein_Quant"

# Setting up the features to use in FeatureMap
ffact = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
fmParams = {}
for k in ffact.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams
exclude = []

fmaps = None

def filter_feature(f):
    if f.GetFamily() in exclude:
        return None
    else:
        return f

def get_raw_features(mol):
    rawFeats = ffact.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filter_feature, rawFeats))
    return filtered

def create_feature_map(mol):
    feats = get_raw_features(mol)
    return FeatMaps.FeatMap(feats=feats, weights=[1]*len(feats),params=fmParams)

def score_featmaps(fm1, fm2):
    """Generate the score for 2 feature maps"""
    if fm1.GetNumFeatures() == 0:
        return 0, 0
    else:
        score = fm1.ScoreFeats(fm2.GetFeatures())
        #utils.log(score, fm1.GetNumFeatures())
        return score, score / fm1.GetNumFeatures()

def get_fmap_scores(mol):
    "Score the molecule using the merged feature map"
    featMap = create_feature_map(mol)
    quantScore, qualScore = score_featmaps(fmaps, featMap)
    return quantScore, qualScore

def build_feat_data(mols):
    "Build the feature maps for the molecules and do the all vs. all comparison"
    featuremaps = []
    scores = []
    for mol1 in mols:
        fm1 = create_feature_map(mol1)
        featuremaps.append(fm1)
        row = []
        for mol2 in mols:
            fm2 = create_feature_map(mol2)
            raw_score, norm_score = score_featmaps(fm1, fm2)
            row.append(norm_score)
            #utils.log(len(data), len(row), raw_score, norm_score)
        scores.append(row)
    return featuremaps, scores

def find_closest(scores):
    #utils.log('Find closest for', len(scores), len(scores[0]))
    best_score = -1.0
    for i in range(len(scores)):
        for j in range(len(scores)):
            if i == j:
                continue
            score = scores[i][j]
            #utils.log('Score:', score)
            if score > best_score:
                best_score = score
                best_row = i
                best_col = j
    return best_score, best_row, best_col

def merge_feat_maps(fmaps, scores):
    "Merge the 2 closest feature maps, remove them form the data and replace with the merged feature map"
    best_score, best_row, best_col = find_closest(scores)
    #utils.log(best_score, best_row, best_col)
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

    #utils.log('Initial:', len(fmaps), len(scores), ','.join([str(len(x)) for x in scores]))
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
        score1 = score_featmaps(fmap, merged)[1]
        score2 = score_featmaps(merged, fmap)[1]
        scores[i].append(score1)
        merged_scores.append(score2)

    fmaps.append(merged)
    merged_scores.append(score_featmaps(merged, merged))
    scores.append(merged_scores)

def score_molecules(molecules, writer):
    total = 0
    success = 0
    errors = 0
    for mol in molecules:
        total += 1
        if mol is None:
            errors += 1
            continue
        try:
            quantScore, qualScore = get_fmap_scores(mol)
            # utils.log('Score:', score)
            if total % 1000 == 0:
                utils.log('Processed molecule', total, '...')
            mol.SetDoubleProp(field_FeatureSteinQualityScore, qualScore)
            mol.SetDoubleProp(field_FeatureSteinQuantityScore, quantScore)
            writer.write(mol)
            success += 1
        except:
            utils.log("Error scoring molecule", sys.exc_info()[0])
            traceback.print_exc()
            errors += 1

    return total, success, errors

def create_featuremap(fragments):

    mols = [m for m in fragments if m]
    fmaps, scores = build_feat_data(mols)
    merged_fmaps = fmaps.copy()
    utils.log('Processing', len(fmaps), 'molecules')
    while len(merged_fmaps) > 1:
        merge_feat_maps(merged_fmaps, scores)
    merged_fmap = merged_fmaps[0]
    utils.log('Created merged feature map with', merged_fmap.GetNumFeatures(), 'features')

    return merged_fmap

### start main execution #########################################

def main():

    # Example usage
    # python -m pipelines.xchem.featurestein_generate_and_score -i ../../data/mpro/poses.sdf.gz -f ../../data/mpro/hits-17.sdf.gz -o output_fs

    global fmaps

    parser = argparse.ArgumentParser(description='FeatureStein scoring with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-f', '--fragments', help='Fragments to use to generate the feature map')
    parser.add_argument('-ff', '--fragments-format', help='Fragments format')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')


    args = parser.parse_args()
    utils.log("FeatureStein Args: ", args)

    source = "featurestein_generate_and_score.py"
    datasetMetaProps = {"source":source, "description": "FeatureStein scoring using RDKit " + rdBase.rdkitVersion}

    clsMappings = {}
    fieldMetaProps = []
    clsMappings[field_FeatureSteinQualityScore] = "java.lang.Float"
    clsMappings[field_FeatureSteinQuantityScore] = "java.lang.Float"
    fieldMetaProps.append({"fieldName":field_FeatureSteinQualityScore,   "values": {"source":source, "description":"FeatureStein quality score"},
                           "fieldName":field_FeatureSteinQuantityScore,   "values": {"source":source, "description":"FeatureStein quantity score"}})

    # generate the feature maps
    frags_input, frags_suppl = rdkit_utils.default_open_input(args.fragments, args.fragments_format)

    fmaps = create_featuremap(frags_suppl)
    frags_input.close()

    # read the ligands to be scored
    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                        'featurestein', args.outformat,
                        valueClassMappings=clsMappings,
                        datasetMetaProps=datasetMetaProps,
                        fieldMetaProps=fieldMetaProps,
                        compress=not args.no_gzip)

    # do the scoring
    total, success, errors = score_molecules(inputs_supplr, writer)
    utils.log('Scored', success, 'molecules.', errors, 'errors.')

    inputs_file.close()
    writer.flush()
    writer.close()
    output.close()

    if args.metrics:
        utils.write_metrics(output_base, {'__InputCount__':total, '__OutputCount__':success, '__ErrorCount__':errors, 'RDKitFeatureMap':success})


if __name__ == "__main__":
    main()
