#!/usr/bin/env nextflow

/* Example Nextflow pipline that runs screen.py followed by conformers.py 
*/

params.qsmiles = 'OC(=O)C1=CC=C(NC2=NC3=C(CN=C(C4=CC(Cl)=CC=C34)C3=C(F)C=CC=C3F)C=N2)C=C1'
params.target = '../../../data/Kinase_inhibs.sdf.gz'
params.simmin = 0.7
params.simmax = 1.0
params.descriptor = 'rdkit'
params.metric = 'tanimoto'
params.num = 1
params.attempts = 0

target = file(params.target)

process rdkitScreen {

	input:
    file target
    file 'screen.py' from file('screen.py')

    output:
    stdout screenOutput
    
    
    """
    ./screen.py --qsmiles '$params.qsmiles' --simmin $params.simmin --simmax $params.simmax -d $params.descriptor -m $params.metric -i $target
    """
    
}

process rdkitConformer {

	input:
    stdin screenOutput
    file 'conformers.py' from file('conformers.py')

    output:
    file 'results.sdf.gz' into results 
    
    
    """
    ./conformers.py -if sdf -n $params.num -a $params.attempts -o results
    """
    
}

results.println { "Results: $it" }
