#!/usr/bin/env nextflow

/* Example Nextflow pipline that runs Docking using rDock 
*/

params.ligands = 'data/hivpr_ligprep_100.sdf.gz'
params.protein = 'data/hivpr_rdock.mol2'
params.prmfile = 'data/hivpr_rdock.prm'
params.asfile =  'data/hivpr_rdock.as'
params.chunk = 25
params.num_dockings = 1
params.top = 1
params.score = null


ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile  = file(params.asfile)

/* Splits the input SD file into multiple files of ${params.chunk} records.
* Each file is sent individually to the ligand_parts channel
*/
process sdsplit {

	input:
    file ligands

    output:
    file 'ligands_part*' into ligand_parts mode flatten
    
    """
	${ligands.name.endsWith('.gz') ? 'gunzip -c' : 'cat'} $ligands | sdsplit -$params.chunk -oligands_part
    """
    
}

/* Docks each file from the ligand_parts channel sending each resulting SD file to the results channel
*/
process rdock {

	input:
    file part from ligand_parts
	file protein
	file prmfile
	file asfile

    output:
    file 'docked.sd' into docked 
    
    
    """
	rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o docked > docked_out.log
    """
    
}

process report {

	input:
	file docked
	
	output:
	file 'results.sd' into results

	"""
	sdsort -n -s -fSCORE docked*.sd | ${params.score == null ? '' : "sdfilter -f'\$SCORE < $params.score' |"} sdfilter -f'\$_COUNT == ${params.top}' > results.sd
	"""

}

/*
 * Collect all hits to a single file called  'all_results'
 */ 
all_results = results.collectFile(name:'all_results.sd')

all_results
    .subscribe { println "Results file: ${it.name}" }
