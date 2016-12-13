#!/usr/bin/env nextflow

/* Example Nextflow pipline that runs screen.py followed by conformers.py 
*/

params.smiles = 'OC(=O)C1=CC=C(NC2=NC3=C(CN=C(C4=CC(Cl)=CC=C34)C3=C(F)C=CC=C3F)C=N2)C=C1'
params.target = '../data/Kinase_inhibs.sdf.gz'
params.simmin = 0.7
params.simmax = 1.0
params.descriptor = 'rdkit'
params.metric = 'tanimoto'
params.num = 1
params.attempts = 0

target = file(params.target)
screenScript = file('screen.py')
conformersScript = file('conformers.py')

process rdkitScreen {

	input:
    file target
    file screenScript

    output:
    stdout screenOutput
    
    
    """
    python $screenScript -smiles '$params.smiles' -simmin $params.simmin -simmax $params.simmax -d $params.descriptor -m $params.metric -i $target
    """
    
}

process rdkitConformer {

	input:
    stdin screenOutput
    file conformersScript

    output:
    file 'results.sdf.gz' into results 
    
    
    """
    python $conformersScript  -n $params.num -a $params.attempts -o results
    """
    
}

results.println { "Results: $it" }
