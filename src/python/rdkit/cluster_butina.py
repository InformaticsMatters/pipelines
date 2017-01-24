#!/usr/bin/env python

import argparse
import utils
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import Butina


descriptors = {
    #'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs':       lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,2,1024),
    'morgan3':     lambda m: AllChem.GetMorganFingerprintAsBitVect(m,3,1024),
    'rdkit':       lambda m: FingerprintMols.FingerprintMol(m),
    #'topo':        lambda m: Torsions.GetTopologicalTorsionFingerprint(m)
}

metrics = {
    'asymmetric':DataStructs.AsymmetricSimilarity,
    'braunblanquet':DataStructs.BulkBraunBlanquetSimilarity,
    'cosine':DataStructs.BulkCosineSimilarity,
    'dice': DataStructs.BulkDiceSimilarity,
    'kulczynski':DataStructs.BulkKulczynskiSimilarity,
    'mcconnaughey':DataStructs.BulkMcConnaugheySimilarity,
    #'onbit':DataStructs.OnBitSimilarity,
    'rogotgoldberg':DataStructs.BulkRogotGoldbergSimilarity,
    'russel':DataStructs.BulkRusselSimilarity,
    'sokal':DataStructs.BulkSokalSimilarity,
    'tanimoto': DataStructs.BulkTanimotoSimilarity
    #'tversky': DataStructs.TverskySimilarity
    }

### start field name defintions #########################################

field_Cluster = "Cluster"

### functions #########################################

def ClusterFps(fps, metric, cutoff):

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        
        func = metrics[metric]
        sims = func(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

def ClustersToMap(clusters):
    d = {}
    i = 0
    for c in clusters:
        for id in c:
            d[id] = i
        i += 1
    return d
    
### start main execution #########################################

def main():

           
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('-t', '--threshold', type=float, default=0.7, help='similarity clustering threshold (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit', help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto', help='similarity metric (default tanimoto)')
    utils.add_default_io_args(parser)
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()
    utils.log("Cluster Args: ",args)

    descriptor = descriptors[args.descriptor]
    if descriptor is None:
        raise ValueError('Invalid descriptor name ' + args.descriptor)

    input,output,suppl,writer,output_base = utils.default_open_input_output(args.input, args.informat, args.output, 'cluster_butina', args.outformat)

    # generate fingerprints
    mols = [x for x in suppl if x is not None]
    fps = [descriptor(x) for x in mols]
    input.close()

    # do clustering
    utils.log("Clustering with descriptor",args.descriptor,"metric",args.metric,"and threshold",args.threshold)
    clusters=ClusterFps(fps, args.metric, 1.0 - args.threshold)
    utils.log("Found",len(clusters),"clusters")
    lookup = ClustersToMap(clusters)

    # write the results
    i = 0
    for mol in mols:
        cluster = lookup[i]
        if (cluster is not None):
            mol.SetIntProp(field_Cluster, cluster)
        i += 1
        writer.write(mol)


    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':i,'RDKitCluster':i})
    
if __name__ == "__main__":
    main()

