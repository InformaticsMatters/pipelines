#!/usr/bin/env python

import utils
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
import gzip, collections
import argparse


### start field name defintions #########################################

field_StructureNum = 'StructureNum'
field_ConformerNum = 'ConformerNum'
field_RMSToCentroid = 'RMSToCentroid'
field_ClusterNum = 'ClusterNum'
field_ClusterCentroid = 'ClusterCentroid'
field_EnergyAbs = 'EnergyAbs'
field_EnergyDelta = 'EnergyDelta'
field_MinimizationConverged = 'MinimizationConverged'

### start function defintions #########################################

def process_mol_conformers(mol, i, numConfs, maxAttempts, pruneRmsThresh, clusterMethod, clusterThreshold, minimizeIterations):
    #utils.log("generating conformers for molecule",i)
    # generate the confomers
    conformerIds = gen_conformers(mol, numConfs, maxAttempts, pruneRmsThresh, True, True, True)
    conformerPropsDict = {}
    minEnergy = 9999999999999
    for conformerId in conformerIds:
        #utils.log("Processing conf",i,conformerId)
        # energy minimise (optional) and energy calculation
        props = collections.OrderedDict()
        energy = calc_energy(mol, conformerId, minimizeIterations, props)
        if energy < minEnergy:
            minEnergy = energy
        conformerPropsDict[conformerId] = props
    # cluster the conformers
    if clusterMethod:
        rmsClusters = cluster_conformers(mol, clusterMethod, clusterThreshold)
        utils.log("Molecule", i, "generated", len(conformerIds), "conformers and", len(rmsClusters), "clusters")
        rmsClustersPerCluster = []
        clusterNumber = 0
    
        for cluster in rmsClusters:
            clusterNumber = clusterNumber+1
            rmsWithinCluster = align_conformers(mol, cluster)
            for conformerId in cluster:
                props = conformerPropsDict[conformerId]
                props[field_ClusterNum] = clusterNumber
                props[field_ClusterCentroid] = cluster[0] + 1
                idx = cluster.index(conformerId)
                if idx > 0:
                    props[field_RMSToCentroid] = rmsWithinCluster[idx-1]
                else:
                    props[field_RMSToCentroid] = 0.0
    else:
        utils.log("Molecule", i, "generated", len(conformerIds), "conformers")
                
    return conformerPropsDict, minEnergy

def write_conformers(mol, i, conformerPropsDict, minEnergy, writer):
    
    for id in range(mol.GetNumConformers()):
        #utils.log("Writing",i,id)
        for name in mol.GetPropNames():
            mol.ClearProp(name)
        mol.SetIntProp(field_StructureNum, i+1)
        mol.SetIntProp(field_ConformerNum, id+1)
        props = conformerPropsDict[id]
        for key in props:
            mol.SetProp(key, str(props[key]))
            energy = props[field_EnergyAbs]
        if energy:
            mol.SetDoubleProp(field_EnergyAbs, energy)
            mol.SetDoubleProp(field_EnergyDelta, energy - minEnergy)
        writer.write(mol, confId=id)

    
def gen_conformers(mol, numConfs=1, maxAttempts=1, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=0)
    #utils.log("generated",len(ids),"conformers")
    return list(ids)
    
def calc_energy(mol, conformerId, minimizeIts, props):
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
    ff.Initialize()
    if minimizeIts > 0:
        props[field_MinimizationConverged] = ff.Minimize(maxIts=minimizeIts)
    e = ff.CalcEnergy()
    props[field_EnergyAbs] = e
    return e
    
def cluster_conformers(mol, mode="RMSD", threshold=2.0):
    if mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(mol)
    else:
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
    #utils.log("generated",len(rms_clusters),"clusters")
    return rms_clusters
    
def align_conformers(mol, clust_ids):
    rmslist = []
    AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
    return rmslist
    
### start main execution #########################################

def main():

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit conformers')
    parser.add_argument('-n', '--num', type=int, default=1, help='number of conformers to generate')
    parser.add_argument('-a', '--attempts', type=int, default=0, help='number of attempts')
    parser.add_argument('-r', '--rmsd', type=float, default=1.0, help='prune RMSD threshold')
    parser.add_argument('-c', '--cluster', type=str.lower, choices=['rmsd', 'tdf'], help='Cluster method (RMSD or TFD). If None then no clustering')
    parser.add_argument('-t', '--threshold', type=float, help='cluster threshold (default of 2.0 for RMSD and 0.3 for TFD)')
    parser.add_argument('-e', '--emin', type=int, default=0, help='energy minimisation iterations (default of 0 means none)')
    utils.add_default_io_args(parser)
    parser.add_argument('--smiles', help='input structure as smiles (incompatible with using files or stdin for input)')

    args = parser.parse_args()

    if not args.threshold:
        if args.cluster == 'tfd':
            args.threshold = 0.3
        else:
            args.threshold = 2.0
        
    utils.log("Conformers Args: ",args)

    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        suppl = [mol]
        input = None
        output,writer,output_base = utils.default_open_output(args.output, 'conformers', args.outformat)
    else:
        input,output,suppl,writer,output_base = utils.default_open_input_output(args.input, args.informat, args.output, 'conformers', args.outformat)

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    if args.meta:
        t = open(output_base + '_types.txt', 'w')
        t.write(field_StructureNum + '=integer\n')
        t.write(field_StructureNum + '=integer\n')
        t.write(field_ConformerNum + '=integer\n')
        t.write(field_EnergyAbs + '=double\n')
        t.write(field_EnergyDelta + '=double\n')
        if args.emin > 0:
            t.write(field_MinimizationConverged + '=boolean\n')
        if args.cluster:
            t.write(field_RMSToCentroid + '=double\n')
            t.write(field_ClusterNum + '=integer\n')
            t.write(field_ClusterCentroid + '=integer\n')
        t.flush()
        t.close()

    i=0
    count=0
    for mol in suppl:
        if mol is None: continue
        m = Chem.AddHs(mol)
        conformerPropsDict, minEnergy = process_mol_conformers(m, i, args.num, args.attempts, args.rmsd, args.cluster, args.threshold, args.emin)
        m = Chem.RemoveHs(m)
        write_conformers(m, i, conformerPropsDict, minEnergy, writer)
        count = count+ m.GetNumConformers()
        i +=1

    if input:
        input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count,'RDKitConformer':count})


if __name__ == "__main__":
    main()

