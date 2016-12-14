#!/usr/bin/env nextflow

params.smiles = 'OC(=O)C1=CC=C(NC2=NC3=C(CN=C(C4=CC(Cl)=CC=C34)C3=C(F)C=CC=C3F)C=N2)C=C1'
params.target = '../data/Kinase_inhibs.sdf.gz'
params.simmin = 0.7
params.simmax = 1.0
params.descriptor = 'rdkit'
params.metric = 'tanimoto'

target = file(params.target)

process rdkitScreen {

	input:
    file target
    file 'screen.py' from file('screen.py')

    output:
    file 'results.sdf.gz' into results 
    
    
    """
    python screen.py -smiles '$params.smiles' -simmin $params.simmin -simmax $params.simmax -d $params.descriptor -m $params.metric -i $target -o results
    """
    
}

results.println { "Results: $it" }
