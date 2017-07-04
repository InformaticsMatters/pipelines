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

import argparse
import gzip
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

from pipelines.utils import utils


### start function definitions #########################################

def fragment(mol, mode, quiet=False):
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol
    else:
        # TODO - handle ties
        if mode == 'hac':
            biggest_count = 0
            for frag in frags:
                hac = frag.GetNumHeavyAtoms()
                if hac > biggest_count:
                    biggest_count = hac
                    biggest_mol = frag
            if not quiet:
                utils.log("Chose fragment from ", len(frags), "based on HAC")
        elif mode == 'mw':
              biggest_mw = 0
              for frag in frags:
                  mw = Descriptors.MolWt(frag)
                  if mw > biggest_mw:
                      biggest_mw = mw
                      biggest_mol = frag
              if not quiet:
                  utils.log("Chose fragment from ", len(frags), "based on MW")
        else:
            raise ValueError('Invalid fragment mode:',mode)
        
        # copy the properties across
        for name in mol.GetPropNames():
            biggest_mol.SetProp(name, mol.GetProp(name))
        return biggest_mol
            
def filter_by_heavy_atom_count(mol, minCount, maxCount, quiet=False):
    hac = mol.GetNumHeavyAtoms()
    if minCount is not None and hac < minCount:
        if not quiet:
            utils.log("HAC", hac, "<", minCount)
        return False
    if maxCount is not None and hac > maxCount:
        if not quiet:
            utils.log("HAC", hac, ">", maxCount)
        return False
    return True

def filter_by_molwt(mol, minMw, maxMw, quiet=False):
    mw = Descriptors.MolWt(mol)
    if minMw is not None and mw < minMw:
        if not quiet:
            utils.log("MolWt", mw, "<", minMw)
        return False
    if maxMw is not None and mw > maxMw:
        if not quiet:
            utils.log("MolWt", mw, ">", maxMw)
        return False
    return True
    
def filter(mol, minHac=None, maxHac=None, minMw=None, maxMw=None, quiet=False):
    if minHac or maxHac:
        if not filter_by_heavy_atom_count(mol, minHac, maxHac, quiet):
            return False
    if minMw or maxMw:
        if not filter_by_molwt(mol, minMw, maxMw, quiet):
            return False 
    return True

### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit filter')
    parser.add_argument('-f', '--fragment', choices=['hac', 'mw'], help='Find single fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight )')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    parser.add_argument('-l', '--limit', type=int, help='Limit output to this many records')
    parser.add_argument('-c', '--chunksize', type=int, help='Split output into chunks of size c. Output will always be files. Names like filter1.sdf.gz, filter2.sdf.gz ...')
    parser.add_argument('-d', '--digits', type=int, default=0, help='When splitting zero pad the file name to this many digits so that they are in sorted order. Names like filter001.sdf.gz, filter002.sdf.gz ...')
    # WARNING: thin output is not appropriate when using --fragment
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode - suppress reporting reason for filtering')
    utils.add_default_io_args(parser)
    args = parser.parse_args()
    utils.log("Filter Args: ", args)
        
    input,suppl = utils.default_open_input(args.input, args.informat)

    
    if args.chunksize:
        chunkNum = 1
        if args.output:
            output_base = args.output
        else:
            output_base = 'filter'
        output_base_chunk = output_base + str(chunkNum).zfill(args.digits)
        output,writer,output_base_chunk = utils.default_open_output(output_base_chunk, output_base_chunk, args.outformat, compress=True)
    else:
        output,writer,output_base_chunk = utils.default_open_output(args.output, "filter", args.outformat, compress=True)
        output_base = output_base_chunk

    utils.log("Writing to " + output_base_chunk)

    i=0
    count = 0
    chunkNum = 1
    for mol in suppl:
        if args.limit and count >= args.limit:
            break
        i +=1
        if mol is None: continue
        if args.fragment:
            mol = fragment(mol, args.fragment, quiet=args.quiet)
        if not filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax, quiet=args.quiet):
            continue
        if args.chunksize:
            if count > 0 and count % args.chunksize == 0:
                writer.close()
                output.close()
                chunkNum += 1
                output_chunk_base = output_base + str(chunkNum).zfill(args.digits)
                utils.log("Writing to " + output_chunk_base)
                output,writer,output_chunk_base = utils.default_open_output(output_chunk_base, output_chunk_base, args.outformat, compress=True)

        count += 1
        writer.write(mol)

    utils.log("Filtered", i, "down to", count, "molecules")
    if args.chunksize:
        utils.log("Wrote", chunkNum, "chunks")
        if (args.digits > 0 and len(str(chunkNum)) > args.digits):
            utils.log("WARNING: not enough digits specified for the number of chunks")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count, 'RDKitFilter':i})

    
if __name__ == "__main__":
    main()

