#!/usr/bin/env nextflow

/* Example Nextflow pipeline that runs rDock docking
*/


params.ligands = 'ligands.sdf.gz'
params.prmfile = 'receptor.prm'
params.protein = 'receptor.mol2'
params.asfile = 'receptor.as'
params.chunk = 25

params.num_dockings = 100
params.top = 1
params.score = null
params.nscore = null
params.limit = 0
params.digits = 4


ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)

/* Splits the input SD file into multiple files of ${params.chunk} records.
* Each file is sent individually to the ligand_parts channel.
* Parts are named so as to be in correct sorted area.
*/
process sdsplit {

    container 'informaticsmatters/rdkit_pipelines:latest'

	input:
    file ligands

    output:
    file 'ligands_part*.sdf' into ligand_parts mode flatten
    
    
    """
    python -m pipelines_utils_rdkit.filter -i $ligands -c $params.chunk -d $params.digits -o ligands_part -of sdf --no-gzip
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
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('ligands', 'docked')[0..-5]} > docked_out.log
    """
}

/* Filter, combine and publish the results
*/
process results {

	input:
	file part from docked_parts.collect()

	output:
	file 'results.sdf' into results

	"""
	sdsort -n -s -fSCORE docked_part*.sd |${params.score == null ? '' : " sdfilter -f'\$SCORE <= $params.score' |"}${params.nscore == null ? '' : " sdfilter -f'\$SCORE.norm <= $params.nscore' |"} sdfilter -f'\$_COUNT <= ${params.top}' > results.sdf
	"""
}
    
results.println { "Results: $it" }
