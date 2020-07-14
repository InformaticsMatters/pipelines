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
See featurestein-generate.py for how to generate the merged feature maps.
"""

from __future__ import print_function
import argparse, os, sys, gzip, pickle, traceback
from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


### start function definitions #########################################

field_FeatureSteinScore = "FeatureStein_Score"

# Setting up the features to use in FeatureMap
ffact = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
fmaps = None


def filter_feature(f):
    result = f.GetFamily() in fmaps.params.keys()
    return result

def get_raw_features(mol):
    rawFeats = ffact.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filter_feature, rawFeats))
    return filtered

def create_feature_map(mol):
    feats = get_raw_features(mol)
    return FeatMaps.FeatMap(feats=feats, weights=[1]*len(feats),params=fmaps.params)

def score_featmaps(fm1):
    "Generate the score for 2 feature maps"
    if fm1.GetNumFeatures() == 0:
        return 0
    else:
        score = fm1.ScoreFeats(fmaps.GetFeatures())
        #utils.log(score, fm1.GetNumFeatures())
        return score / fm1.GetNumFeatures()

def get_fmap_score(mol):
    featMap = create_feature_map(mol)
    score = score_featmaps(featMap)
    return score

def process(inputs, writer):
    total = 0
    success = 0
    errors = 0
    for mol in inputs:
        total += 1
        if mol is None:
            errors += 1
            continue
        try:
            score = get_fmap_score(mol)
            # utils.log('Score:', score)
            if total % 1000 == 0:
                utils.log('Processed molecule', total, '...')
            mol.SetDoubleProp(field_FeatureSteinScore, score)
            writer.write(mol)
            success += 1
        except:
            utils.log("Error scoring molecule", sys.exc_info()[0])
            traceback.print_exc()
            errors += 1

    return total, success, errors

### start main execution #########################################

def main():

    # Example usage
    # python -m pipelines.xchem.featurestein-score -i ../../data/mpro/poses.sdf.gz -f mpro-fstein.p -o fstein

    global fmaps

    parser = argparse.ArgumentParser(description='FeatureStein scoring with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-f', '--feat-map', help='Feature Map pickle to score with')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')


    args = parser.parse_args()
    utils.log("FeatureStein Args: ", args)

    source = "featurestein-score.py"
    datasetMetaProps = {"source":source, "description": "FeatureStein scoring using RDKit " + rdBase.rdkitVersion}

    clsMappings = {}
    fieldMetaProps = []
    clsMappings[field_FeatureSteinScore] = "java.lang.Float"
    fieldMetaProps.append({"fieldName":field_FeatureSteinScore,   "values": {"source":source, "description":"FeatureStein score"}})

    pkl_file = open(args.feat_map, 'rb')
    fmaps = pickle.load(pkl_file)
    utils.log('FeatureMap has', fmaps.GetNumFeatures(), "features")

    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                        'featurestein', args.outformat,
                        valueClassMappings=clsMappings,
                        datasetMetaProps=datasetMetaProps,
                        fieldMetaProps=fieldMetaProps,
                        compress=not args.no_gzip)

    # this does the processing
    total, success, errors = process(inputs_supplr, writer)

    inputs_file.close()
    writer.flush()
    writer.close()
    output.close()

    if args.metrics:
        utils.write_metrics(output_base, {'__InputCount__':total, '__OutputCount__':success, '__ErrorCount__':errors, 'RDKitFeatureMap':success})


if __name__ == "__main__":
    main()
