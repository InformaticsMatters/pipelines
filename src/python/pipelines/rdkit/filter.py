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
from rdkit.Chem import Descriptors
from pipelines.utils import utils
from pipelines.rdkit import mol_utils


### start function definitions #########################################


            
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
    parser.add_argument('-f', '--fragment', choices=['hac', 'mw'], help='Find single fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    parser.add_argument('-l', '--limit', type=int, help='Limit output to this many records')
    parser.add_argument('-c', '--chunksize', type=int, help='Split output into chunks of size c. Output will always be files. Names like filter1.sdf.gz, filter2.sdf.gz ...')
    parser.add_argument('-d', '--digits', type=int, default=0, help='When splitting zero pad the file name to this many digits so that they are in sorted order. Names like filter001.sdf.gz, filter002.sdf.gz ...')
    parser.add_argument('-r', '--rename', action='append', help='Rename field (fromname:toname)')
    parser.add_argument('--delete', action='append', help='Delete field')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    # WARNING: thin output is not appropriate when using --fragment
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode - suppress reporting reason for filtering')
    utils.add_default_io_args(parser)
    args = parser.parse_args()
    utils.log("Filter Args: ", args)

    field_renames = {}
    if args.rename:
        for t in args.rename:
            parts = t.split(':')
            if len(parts) != 2:
                raise ValueError('Invalid field rename argument:',t)
            field_renames[parts[0]] = parts[1]
    if args.delete:
        for f in args.delete:
            field_renames[f] = None
        
    input,suppl = utils.default_open_input(args.input, args.informat)

    if args.chunksize:
        chunkNum = 1
        if args.output:
            output_base = args.output
        else:
            output_base = 'filter'
        output_base_chunk = output_base + str(chunkNum).zfill(args.digits)
        output,writer,output_base_chunk = utils.default_open_output(output_base_chunk, output_base_chunk, args.outformat, compress=not args.no_gzip)
    else:
        output,writer,output_base_chunk = utils.default_open_output(args.output, "filter", args.outformat, compress=not args.no_gzip)
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
            mol = mol_utils.fragment(mol, args.fragment, quiet=args.quiet)
        if not filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax, quiet=args.quiet):
            continue
        if args.chunksize:
            if count > 0 and count % args.chunksize == 0:
                # new chunk, so create new writer
                writer.close()
                output.close()
                chunkNum += 1
                output_chunk_base = output_base + str(chunkNum).zfill(args.digits)
                utils.log("Writing to " + output_chunk_base)
                output,writer,output_chunk_base = utils.default_open_output(output_chunk_base, output_chunk_base, args.outformat, compress=not args.no_gzip)

        for from_name in field_renames:
            to_name = field_renames[from_name]
            if mol.HasProp(from_name):
                val = mol.GetProp(from_name)
                mol.ClearProp(from_name)
                if to_name:
                    mol.SetProp(to_name, val)

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

