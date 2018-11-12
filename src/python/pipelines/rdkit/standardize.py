#!/usr/bin/env python

# Copyright 2018 Informatics Matters Ltd.
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

import argparse

from rdkit import DataStructs, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils


### functions #########################################

#lfc = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()


def standardize(mol, neutralise, fragment):
    """

    :param mol: The molecule to standardize
    :param neutralise: Boolean for whether to neutralise the molecule
    :param fragment: The approach for choosing the largest fragment. Either 'hac' or 'mw'. If not specified the whole
    molecule is used.
    :return: The standardized molecule
    """
    mol = rdMolStandardize.Cleanup(mol)
    #mol = lfc.choose(mol)
    # We use our own largest fragment picker as the RDKit one behaves slightly differently
    if fragment:
        mol = mol_utils.fragment(mol, fragment)
    if neutralise:
        mol = uncharger.uncharge(mol)
    return mol


### start main execution #########################################

def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='RDKit Standardize')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'], help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--neutralise', action='store_true', help='Neutralise the molecule')

    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()
    utils.log("Standardize Args: ", args)

    # handle metadata
    source = "standardize.py"
    datasetMetaProps = {"source":source, "description": "Standardize using RDKit " + rdBase.rdkitVersion}
    clsMappings = {}
    fieldMetaProps = []


    input,output,suppl,writer,output_base = rdkit_utils.\
        default_open_input_output(args.input, args.informat, args.output,
                                  'standardize', args.outformat,
                                  thinOutput=False, valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)
    i = 0
    total = 0
    for mol in suppl:
        if mol is None:
            i += 1
            continue
        m = standardize(mol, args.neutralise, args.fragment_method)
        writer.write(m)
        total += 1

    input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':total, 'RDKitStandardize':i})

if __name__ == "__main__":
    main()

