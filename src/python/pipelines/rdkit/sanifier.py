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

### Use MolVS to do tautomer enumeration, sterochemistry enumeration, charge neutralisation.

from pipelines.utils import utils
import argparse
from sanify_utils import enumerateStereoIsomers,enumerateTautomers,getCanonTautomer,STANDARD_MOL_METHODS


def write_out(mols,count,writer):
    for mol in mols:
        count += 1
        if mol is None: continue
        writer.write(mol)
    return count

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit molecule standardiser / enumerator')
    utils.add_default_io_args(parser)
    parser.add_argument('-et', '--enumerate_tauts', default=False, help='Enumarate all tautomers', type=bool)
    parser.add_argument('-es', '--enumerate_stereo', default=False, help='Enumarate all stereochem', type=bool)
    parser.add_argument('-st', '--standardize', default=False, help='Standardize molecules. Cannot  be true if enumerate is on.', type=bool)
    parser.add_argument('-stm','--standardize_method', default="molvs",choices=STANDARD_MOL_METHODS.keys(),help="Chose the method to standardize.")
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()

    if args.standardize and args.enumerate_tauts:
        raise ValueError("Cannot Enumerate Tautomers and Standardize")

    if args.standardize and args.enumerate_stereo:
        raise ValueError("Cannot Enumerate Stereo and Standardize")

    if args.standardize:
        getStandardMolecule = STANDARD_MOL_METHODS[args.standardize_method]

    input ,output ,suppl ,writer ,output_base = utils.default_open_input_output(args.input, args.informat, args.output, 'screen', args.outformat, thinOutput=args.thin)
    i=0
    count = 0
    for mol in suppl:
        i +=1
        if mol is None: continue

        if args.standardize:
            count = write_out([getStandardMolecule(mol)],count,writer)
        if args.enumerate_stereo:
            count = write_out(enumerateStereoIsomers(mol),count,writer)
        if args.enumerate_tauts:
            count = write_out(enumerateTautomers(mol),count,writer)


    utils.log("Cleaned up "+str(i)+" molecules, resulting in "+str(count)+" outputs")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count , 'RDKitScreen':count })

    return count

if __name__ == "__main__":
    main()
