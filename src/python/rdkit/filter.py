#!/usr/bin/env python

import utils
import sys, gzip
from rdkit import Chem
from rdkit.Chem import Descriptors
import argparse

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
            biggest_mol.PutProp(name, mol.GetProp(name))
        return biggest_mol
            
def filter_by_heavy_atom_count(mol, minCount, maxCount, quiet=False):
    hac = mol.GetNumHeavyAtoms()
    if minCount is not None and hac < minCount:
        if not quiet:
            utils.log("HAC",hac,"<", minCount)
        return False
    if maxCount is not None and hac > maxCount:
        if not quiet:
            utils.log("HAC",hac,">", maxCount)
        return False
    return True

def filter_by_molwt(mol, minMw, maxMw, quiet=False):
    mw = Descriptors.MolWt(mol)
    if minMw is not None and mw < minMw:
        if not quiet:
            utils.log("MolWt",mw,"<",minMw)
        return False
    if maxMw is not None and mw > maxMw:
        if not quiet:
            utils.log("MolWt",mw,">", maxMw)
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
    parser.add_argument('-c', '--chunksize', type=int, help='Split output into chunks of size c. Output will always be files. Names like filter01.sdf.gz, filter02.sdf.gz ...')
    # WARNING: thin output is not appropriate when using --fragment
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode - suppress reporting reason for filtering')
    utils.add_default_io_args(parser)
    args = parser.parse_args()
    utils.log("Filter Args: ",args)
        
    input,suppl = utils.default_open_input(args.input, args.informat)
    if args.output:
        output_base = args.output
    else:
        output_base = 'filter'

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    #t = open(output_base + '_types.txt', 'w')
    #t.write(field_Similarity + '=integer\n')
    #t.flush()
    #t.close()
    
    if args.chunksize:
        chunkNum = 1
        output = gzip.open(output_base + str(chunkNum) + '.sdf.gz','w+')
    elif args.output:
        output = gzip.open(output_base + '.sdf.gz','w+')
    else:
        output = sys.stdout
        
    writer = Chem.SDWriter(output)

    i=0
    count = 0
    chunkCount = 1
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
                chunkCount += 1
                output = gzip.open(output_base + str(chunkCount) + '.sdf.gz','w+')
                writer = Chem.SDWriter(output)

        count += 1
        writer.write(mol)

    utils.log("Filtered",i,"down to",count,"molecules")
    if args.chunksize:
        utils.log("Wrote",chunkCount,"chunks")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count,'RDKitFilter':i})

    
if __name__ == "__main__":
    main()

