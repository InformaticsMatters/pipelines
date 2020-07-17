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
Expand a set of molecules using the fragment network
"""

from rdkit import Chem, rdBase

import os, argparse
import requests

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils


def process(inputs, hac_min=None, hac_max=None, rac_min=None, rac_max=None, hops=1,
            server='https://fragnet-external.xchem-dev.diamond.ac.uk', token=None,
            index_as_filename=False):

    count = 0
    for mol in inputs:
        count += 1
        name = mol.GetProp('_Name')
        smiles = Chem.MolToSmiles(mol)
        utils.log('Processing', name, smiles)
        url = server + '/fragnet-search/rest/v2/search/expand/' + smiles
        params = {'hops': hops, 'hac_min': hac_min, 'hac_max': hac_max, 'rac_min': rac_min, 'rac_max': rac_max}

        headers = {}
        if token:
            headers['Authorization'] = 'bearer ' + token

        r = requests.get(url, params=params, headers=headers)
        if r.status_code == requests.codes.ok:
            j = r.json()
            num_mols = j['size']
            utils.log('JSON OK', num_mols, 'results')
            basename = None
            if not index_as_filename:
                if mol.HasProp('_Name'):
                    name = mol.GetProp('_Name')
                    if name:
                        basename = name
            if not basename:
                basename = str(count)

            Chem.MolToMolFile(mol, basename + '.mol')

            with open(basename + '.smi', 'w') as f:
                members = j['members']
                for member in members:
                    f.write(member['smiles'] + '\n')
        else:
            utils.log('Request failed with status code', r.status_code)


def main():
    # Example usage:
    # 1. Create keycloak token:
    # export KEYCLOAK_TOKEN=$(curl -d "grant_type=password" -d "client_id=fragnet-search" -d "username=<username>" -d "password=<password>" \
    #   https://squonk.it/auth/realms/squonk/protocol/openid-connect/token 2> /dev/null | jq -r '.access_token')
    #
    # 2. Run the module:
    #  python -m pipelines.xchem.fragnet_expand -i ../../data/mpro/hits-17.sdf.gz --token $KEYCLOAK_TOKEN

    parser = argparse.ArgumentParser(description='Fragnet expand scoring with RDKit')
    parameter_utils.add_default_input_args(parser)
    parser.add_argument('--hac-min', type=int, default=3, help='The min change in heavy atom count')
    parser.add_argument('--hac-max', type=int, default=3, help='The max change in heavy atom count')
    parser.add_argument('--rac-min', type=int, default=1, help='The min change in ring atom count')
    parser.add_argument('--rac-max', type=int, default=1, help='The max change in ring atom count')
    parser.add_argument('--hops', type=int, default=1, help='The number of graph traversals (hops)')
    parser.add_argument('-s', '--server', default='https://fragnet-external.xchem-dev.diamond.ac.uk', help='The fragnet search server')
    parser.add_argument('--token', help='Keycloak auth token (or specify as KEYCLOAK_TOKEN env variable')
    parser.add_argument('--index-as-filename', action='store_true', help='Use the index as the file name instead of the molecule name')


    args = parser.parse_args()
    utils.log("FragnetExpand Args: ", args)

    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)

    if args.token:
        auth_token = args.token
    else:
        auth_token = os.getenv('KEYCLOAK_TOKEN')
    if not auth_token:
        utils.log('WARNING: not authentication token found in environment variable KEYCLOAK_TOKEN')

    # this does the processing
    process(inputs_supplr, hac_min=args.hac_min, hac_max=args.hac_max, rac_min=args.rac_min, rac_max=args.rac_max,
            hops=args.hops,
            server=args.server, token=auth_token,
            index_as_filename=args.index_as_filename)

    inputs_file.close()


if __name__ == "__main__":
    main()
