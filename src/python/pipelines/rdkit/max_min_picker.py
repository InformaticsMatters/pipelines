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

import argparse, logging, time

from rdkit import DataStructs, rdBase, SimDivFilters
from rdkit.Chem import AllChem, MACCSkeys

from pipelines.utils import utils
from pipelines.rdkit import mol_utils

descriptors = {
    #'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,2,1024),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3,1024)
}


### functions #########################################


def performPick(fpBitVector, howManyToPick, similarityThreshold, firstPicks):
    picker = SimDivFilters.MaxMinPicker()
    # LazyPickWithThreshold( (MaxMinPicker)self, (AtomPairsParameters)distFunc, (int)poolSize, (int)pickSize, (float)threshold [, (AtomPairsParameters)firstPicks=() [, (int)seed=-1]]) -> tuple :
    picks, thresh = picker.LazyBitVectorPickWithThreshold(fpBitVector, len(fpBitVector), howManyToPick, similarityThreshold, firstPicks=firstPicks)
    return picks, thresh


### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit Butina Cluster')
    parser.add_argument('-t', '--threshold', type=float, default=0.0, help='similarity threshold (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='morgan2', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-n', '--num', type=int, help='maximum number to pick for diverse subset selection')
    parser.add_argument('-s', '--seed-molecules', help='optional file containing any seed molecules that have already been picked')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], default='hac', help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--output-fragment', action='store_true', help='Output the biggest fragment rather than the original molecule')
    utils.add_default_io_args(parser)

    args = parser.parse_args()
    utils.log("MaxMinPicker Args: ", args)

    descriptor = descriptors[args.descriptor]
    if descriptor is None:
        raise ValueError('No descriptor specified')

    if not args.num and not args.threshold:
        raise ValueError('--num or --threshold arguments must be specified, or both')

    # handle metadata
    source = "max_min_picker.py"
    datasetMetaProps = {"source":source, "description": "MaxMinPicker using RDKit " + rdBase.rdkitVersion}

    ### generate fingerprints
    fps = []
    mols = []
    errors = 0

    # first the initial seeds, if specified
    firstPicks = []
    num_seeds = 0
    if args.seed_molecules:
        seedsInput,seedsSuppl = utils.default_open_input(args.seed_molecules, None)
        start = time.time()
        errors += mol_utils.fragmentAndFingerprint(seedsSuppl, mols, fps, descriptor, fragmentMethod=args.fragment_method, outputFragment=args.output_fragment, quiet=args.quiet)
        end = time.time()
        seedsInput.close()
        num_seeds = len(fps)
        utils.log("Read", len(fps), "fingerprints for seeds in", end-start, "secs,", errors, "errors")
        firstPicks = list(range(num_seeds))

    # now the molecules to pick from
    input,output,suppl,writer,output_base = utils.default_open_input_output(args.input, args.informat, args.output, 'cluster_butina',
                                                                            args.outformat, datasetMetaProps=datasetMetaProps)
    # reset the mols list as we don't need the seeds, only the candidates
    mols = []
    start = time.time()
    errs = mol_utils.fragmentAndFingerprint(suppl, mols, fps, descriptor, fragmentMethod=args.fragment_method, outputFragment=args.output_fragment, quiet=args.quiet)
    end = time.time()
    errors += errs

    input.close()
    num_fps = len(fps)
    num_candidates = num_fps - num_seeds
    utils.log("Read", num_candidates, "fingerprints for candidates in", end-start, "secs,", errs, "errors")

    if not args.num:
        num_to_pick = num_candidates
    elif args.num > num_candidates:
        num_to_pick = num_candidates
        utils.log("WARNING: --num argument (", args.num, ") is larger than the total number of candidates (", num_candidates, ") - resetting to", num_candidates)
    else:
        num_to_pick = args.num

    ### do picking
    utils.log("MaxMinPicking with descriptor", args.descriptor, "and threshold", args.threshold, ",", num_seeds, "seeds,", num_candidates, "candidates", num_fps, "total")
    start = time.time()
    picks, thresh = performPick(fps, num_to_pick + num_seeds, args.threshold, firstPicks)
    end = time.time()
    num_picks = len(picks)

    utils.log("Found", num_picks, "molecules in", end-start, "secs, final threshold", thresh)
    utils.log("Picks:", list(picks[num_seeds:]))
    del fps

    # we want to return the results in the order they were in the input so first we record the order in the pick list
    indices = {}
    i = 0
    for idx in picks[num_seeds:]:
        indices[idx] = i
        i += 1
    # now do the sort
    sorted_picks = sorted(picks[num_seeds:])
    # now write out the mols in the correct order recording the value in the pick list as the PickIndex property
    i = 0
    for idx in sorted_picks:
        mol = mols[idx - num_seeds] # mols array only contains the candidates
        mol.SetIntProp("PickIndex", indices[idx] + 1)
        writer.write(mol)
        i += 1
    utils.log("Output", i, "molecules")

    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        metrics = {}
        status_str = "{} compounds picked. Final threshold was {}.".format(i, thresh)
        if errors > 0:
            metrics['__ErrorCount__'] = errors
            status_str = status_str + " {} errors.".format(errors)

        metrics['__StatusMessage__'] = status_str
        metrics['__InputCount__'] = num_fps
        metrics['__OutputCount__'] = i
        metrics['RDKitMaxMinPicker'] = num_picks

        utils.write_metrics(output_base, metrics)

if __name__ == "__main__":
    main()
