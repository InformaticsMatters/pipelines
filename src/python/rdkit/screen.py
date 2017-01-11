#!/usr/bin/env python

import utils, filter
import sys, gzip, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

### start field name defintions #########################################

field_Similarity = "Similarity"

### start main execution #########################################

descriptors = {
    #'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2':     lambda m: AllChem.GetMorganFingerprint(m,2),
    'morgan3':     lambda m: AllChem.GetMorganFingerprint(m,3),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m),
    #'topo':        lambda m: Torsions.GetTopologicalTorsionFingerprint(m)
}

metrics = {
    'asymmetric':DataStructs.AsymmetricSimilarity,
    'braunblanquet':DataStructs.BraunBlanquetSimilarity,
    'cosine':DataStructs.CosineSimilarity,
    'dice': DataStructs.DiceSimilarity,
    'kulczynski':DataStructs.KulczynskiSimilarity,
    'mcconnaughey':DataStructs.McConnaugheySimilarity,
    #'onbit':DataStructs.OnBitSimilarity,
    'rogotgoldberg':DataStructs.RogotGoldbergSimilarity,
    'russel':DataStructs.RusselSimilarity,
    'sokal':DataStructs.SokalSimilarity,
    'tanimoto': DataStructs.TanimotoSimilarity
    #'tversky': DataStructs.TverskySimilarity
}

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('--smiles', help='query structure as smiles (incompatible with -molfile arg)')
    parser.add_argument('--molfile', help='query structure as filename in molfile format (incompatible with -smiles arg)')
    parser.add_argument('--simmin', type=float, default=0.7, help='similarity lower cutoff (1.0 means identical)')
    parser.add_argument('--simmax', type=float, default=1.0, help='similarity upper cutoff (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', choices=list(descriptors.keys()), default='rdkit', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric',
                    choices=list(metrics.keys()),
                    default='tanimoto', help='similarity metric (default tanimoto)')
    parser.add_argument('-f', '--fragment', choices=['hac', 'mw'], help='Find single fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight )')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    args = parser.parse_args()
    utils.log("Screen Args: ",args)

    descriptor = descriptors[args.descriptor.lower()]
    metric = metrics[args.metric.lower()]

    
    if args.smiles and args.molfile:
        raise ValueError('Cannot specify -smiles and -molfile arguments together')
    elif args.smiles:
        query_rdkitmol = Chem.MolFromSmiles(args.smiles)
    elif args.molfile:
        query_rdkitmol = Chem.MolFromMolFile(args.molfile)
    else:
        raise ValueError('No query structure specified')
    
    query_fp = descriptor(query_rdkitmol)

    input,output,suppl,writer,output_base = utils.default_open_input_output(args.input, args.informat, args.output, 'screen', args.outformat)

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    if args.meta:
        t = open(output_base + '_types.txt', 'w')
        t.write(field_Similarity + '=integer\n')
        t.flush()
        t.close()

    i=0
    count = 0
    for mol in suppl:
        i +=1
        if mol is None: continue
        if args.fragment:
            mol = filter.fragment(mol, args.fragment, quiet=args.quiet)
        if not filter.filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax, quiet=args.quiet):
            continue
        target_fp = descriptor(mol)
        sim = metric(query_fp, target_fp)
    
        if sim >= args.simmin and sim <= args.simmax:
            count +=1
            if not args.quiet:
                utils.log(i,sim)
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            mol.SetDoubleProp(field_Similarity, sim)
            writer.write(mol)

    utils.log("Found",count,"similar molecules")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count,'RDKitScreen':count})

    return count
    
if __name__ == "__main__":
    main()

