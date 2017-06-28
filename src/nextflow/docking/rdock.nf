#!/usr/bin/env nextflow

/* Example Nextflow pipline that runs Docking using rDock 
*/

params.prmfile = 'receptor.prm'
params.ligands = 'ligands.sdf.gz'
params.protein = 'receptor.mol2'
params.asfile =  'receptor.as'
params.chunk = 25
params.num_dockings = 100
params.top = 1
params.score = null

prmfile = file(params.prmfile)
ligands = file(params.ligands)
protein = file(params.protein)
asfile  = file(params.asfile)

/* Splits the input SD file into multiple files of ${params.chunk} records.
* Each file is sent individually to the ligand_parts channel.
* Parts are renamed so as to be in correct sorted area.
* TODO - the renaming is a bit hacky - look for a better solution 
*/
process sdsplit {

	input:
    file ligands

    output:
    file 'ligands_part*' into ligand_parts mode flatten
    
    
    """ 
	${ligands.name.endsWith('.gz') ? "zcat $ligands > ${ligands.name[0..-4]} && " : ''}sdsplit -$params.chunk -oligands_part ${ligands.name.endsWith('.gz') ? ligands.name[0..-4] : ligands.name} > sdsplit.log
	files=(ligands_part?.sd) && if [ -e "\${files[0]}" ]; then rename ligands_part ligands_part0 ligands_part?.sd; fi
	files=(ligands_part??.sd) && if [ -e "\${files[0]}" ]; then rename ligands_part ligands_part00 ligands_part??.sd; fi
	files=(ligands_part???.sd) && if [ -e "\${files[0]}" ]; then rename ligands_part ligands_part000 ligands_part???.sd; fi
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
    file 'docked_part*.sd' into docked_parts 
    
    """
	rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('ligands', 'docked')[0..-4]} > docked_out.log
    """
}

/* Filter, combine and publish the results
*/
process results {

	publishDir './', mode: 'copy'
	
	input:
	file ligands
	file part from docked_parts.collect()
	
	output:
	file 'output.sdf.gz'
	file 'output_metrics.txt'
	
 
	"""
	sdsort -n -s -fSCORE docked_part*.sd | ${params.score == null ? '' : "sdfilter -f'\$SCORE < $params.score' |"} sdfilter -f'\$_COUNT <= ${params.top}' | gzip > output.sdf.gz
	echo -n 'DockingRDock=' > output_metrics.txt
	echo \$((`${ligands.name.endsWith('.gz') ? 'zcat' : 'cat'} $ligands.name | fgrep -c '\$\$\$\$'` * $params.num_dockings)) >> output_metrics.txt
	echo -n '__OutputCount__=' >> output_metrics.txt
	zcat output.sdf.gz | fgrep -c '\$\$\$\$' >> output_metrics.txt
	"""
}

    
    
