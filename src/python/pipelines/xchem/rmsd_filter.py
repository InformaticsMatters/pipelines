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
Filter poses using a RMSD cuttoff. The poses for an individual molecule are expected to be sequential records, identified
by the field name specified by the --field argument, and are expected to be already ranked best to worst within each group.
Within each group the first (best) molecule is kept and then the remaining molecules are compared to those already accepted
in the group and if the RMSD is less than the threshold specified by the --cutoff-rmsd argument for any of the already
accepted molecules then the molecule is discarded. Only if its RMSD is greater than the threshold for all already accepted
members of the group is the molecule retained.
"""

from __future__ import print_function
import argparse, os, sys, gzip, traceback
from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


### start function definitions #########################################

def process(supplr, writer, grouping_field, rms_cutoff):
    molconfs = None
    mols_to_keep = []
    group = None
    count = 0
    errors = 0
    kept = 0
    confs = 0
    groups = 0
    for mol in supplr:
        count += 1
        confs += 1
        if not mol:
            errors += 1
            continue
        g = mol.GetProp(grouping_field)
        if not molconfs:
            utils.log('New group', g)
            molconfs = mol
            mols_to_keep.append(mol)
            kept += 1
            confs = 1
            groups += 1
            group = mol.GetProp(grouping_field)
        else:
            if group == g:
                conf_id = molconfs.AddConformer(mol.GetConformer(-1), assignId=True)
                # utils.log('  added conf', conf_id)
                num_confs = molconfs.GetNumConformers()
                # utils.log('  Now have', num_confs, 'conformers')
                keep = True
                for i in range(num_confs - 1):
                    rms = AllChem.GetConformerRMS(molconfs, i, num_confs - 1, prealigned=True)
                    # utils.log('  comparing', i, num_confs -1, rms)
                    if rms < rms_cutoff:
                        molconfs.RemoveConformer(conf_id)
                        # utils.log('  dropping conformer', conf_id, 'with RMS', rms)
                        keep = False
                        break
                if keep:
                    mols_to_keep.append(mol)
                    kept += 1
                    # utils.log('  keeping conformer', conf_id, 'now have', len(mols_to_keep))

            else:
                utils.log('  keeping', len(mols_to_keep), 'out of', confs, 'conformers')
                write_mols(writer, mols_to_keep)
                utils.log('New group', g)
                molconfs = mol
                mols_to_keep = [mol]
                kept += 1
                confs = 1
                group = g
                groups += 1

    # handle the final group
    if mols_to_keep:
        utils.log('  keeping', len(mols_to_keep), 'out of', confs, 'conformers')
        write_mols(writer, mols_to_keep)

    return count, groups, kept, errors


def write_mols(writer, mols):
    for mol in mols:
        writer.write(mol)
    writer.flush()


### start main execution #########################################

def main():

    # Example usage
    # python -m pipelines.xchem.rmsd_filter -i ../../data/mpro/poses.sdf.gz -o output -c 0.5

    parser = argparse.ArgumentParser(description='RSMD filter')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-c', '--cutoff-rmsd', type=float, help='RMSD cutoff')
    parser.add_argument('-f', '--field',  default='_Name', help='Field to group records')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--metrics', action='store_true', help='Write metrics')


    args = parser.parse_args()
    utils.log("RSMD filter Args: ", args)

    source = "rmsd_filter.py"
    datasetMetaProps = {"source":source, "description": "RMSD filter " + rdBase.rdkitVersion}

    clsMappings = {}
    fieldMetaProps = []
    # clsMappings[field_FeatureSteinQualityScore] = "java.lang.Float"
    # clsMappings[field_FeatureSteinQuantityScore] = "java.lang.Float"
    # fieldMetaProps.append({"fieldName":field_FeatureSteinQualityScore,   "values": {"source":source, "description":"FeatureStein quality score"},
    #                        "fieldName":field_FeatureSteinQuantityScore,   "values": {"source":source, "description":"FeatureStein quantity score"}})


    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                        'rmsd_filter', args.outformat,
                        valueClassMappings=clsMappings,
                        datasetMetaProps=datasetMetaProps,
                        fieldMetaProps=fieldMetaProps,
                        compress=not args.no_gzip)

    # this does the processing
    count, groups, kept, errors = process(inputs_supplr, writer, args.field, args.cutoff_rmsd)
    utils.log('Processing complete.', count, 'records processed with', groups, 'groups.', kept, 'records retained.', errors, 'errors')

    inputs_file.close()
    writer.flush()
    writer.close()
    output.close()

    if args.metrics:
        utils.write_metrics(output_base, {'__InputCount__':total, '__OutputCount__':success, '__ErrorCount__':errors, 'RDKitFeatureMap':success})


if __name__ == "__main__":
    main()
