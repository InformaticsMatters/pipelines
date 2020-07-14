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

import argparse, os, sys, json, traceback
from pipelines_utils import utils
from pipelines_utils import utils

def gen_filename(id, generate_filenames):
    if generate_filenames:
        return str(count)
    else:
        return id

def execute(candidates_json, generate_filenames):

    with open(candidates_json, 'r') as f:
        candidates = json.load(f)
        queries = candidates['queries']['molecules']
        results = candidates['results']
        hitCounts = candidates['hitCounts']
        utils.log('Processing', len(queries), 'queries and', len(results), 'results')

        num_mols = 0
        num_hits = 0

        count = 0
        ids2Filenames = {}
        for query in queries:
            id = query['id']
            if id in hitCounts:
                molfile = query['originalMol']
                if generate_filenames:
                    fname = str(count).zfil(3)
                else:
                    fname = id
                utils.log('Using file name of', fname)

                with open(fname + '.mol', 'w') as f:
                    f.write(molfile)
                    num_hits += 1
                ids2Filenames[id] = fname
                count += 1

        writers = {}

        for result in results:

            num_mols += 1

            for id in result['sourceMols']:

                if id in writers:
                    writer =  writers[id]
                else:
                    fname = ids2Filenames[id]
                    writer = open(fname + '.smi', 'w')
                    writers[id] = writer

                smiles = result['smiles']
                #utils.log('Processing', smiles)

                writer.write(smiles + '\n')

        for w in writers.values():
            w.close()

        utils.log('Totals - hits:', num_hits, 'outputs:', num_mols)

def main():
    """
    Example usage:
    python -m pipelines.xchem.split-fragnet-candidates -i ../../data/mpro/expanded-17.json

    :return:
    """

    parser = argparse.ArgumentParser(description='Split fragnet candidates - Split fragment network expansion into individual sets')

    parser.add_argument('-i', '--input', help='JSON containing the expanded candidates)')
    parser.add_argument('-g', '--generate-filenames', action='store_true', help='Use automatically generated file names instead of the title field)')

    args = parser.parse_args()
    utils.log("Split fragnet candidates args: ", args)

    infile = args.input

    execute(infile, args.generate_filenames)

if __name__ == "__main__":
    main()
