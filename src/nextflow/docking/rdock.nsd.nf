#!/usr/bin/env nextflow

/* Squonk Nextflow pipline that runs Docking using rDock.
* The contents of the zip file specified by params.receptor must contain the following:
* 1. receptor.mol2 - the prepared protein in mol2 format
* 2. receptor.as - the receptor active site definition
* 3. receptor.prm - the rDock configuration file that refers to receptor.mol2 in its RECEPTOR_FILE property.
* This zip file is unzipped and the contents used by rDock.
* To test this manually run something like this:
* nextflow run src/nextflow/docking/rdock.nsd.nf --ligands data/dhfr_3d.data.gz --receptor data/hivpr.config.zip --num_dockings 5 -with-docker informaticsmatters/rdkit_pipelines
*/

params.ligands = "$baseDir/ligands.data.gz"
params.receptor = "$baseDir/config.zip"
params.chunk = 25
params.num_dockings = 100
params.top = 1
params.score = null
params.nscore = null
params.limit = 0
params.digits = 4


ligands = file(params.ligands)
receptorzip = file(params.receptor)

process unzip_config {

    beforeScript 'chmod g+w .'
    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file receptorzip

    output:
    file 'receptor.prm' into prmfile
    file 'receptor.mol2' into protein
    file 'receptor.as' into asfile

    """
    unzip $receptorzip
    """

}

/* Splits the input into multiple files of ${params.chunk} records.
*/
process splitter {

    beforeScript 'chmod g+w .'
    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file ligands

    output:
    file 'ligands_part*.sdf' into ligands_parts mode flatten
    file 'ligands_part_metrics.txt' into splitter_metrics

    """
    python -m pipelines_utils_rdkit.filter -i $ligands -c $params.chunk -l $params.limit -d $params.digits -o ligands_part -of sdf --no-gzip --meta
    """
}

/* Docks each file from the ligand_parts channel sending each resulting SD file to the results channel
*/
process rdock {

    container 'informaticsmatters/rdock-mini:latest'
    // change permissions on the work dir so that the rdock user in the container
    // can write to the directory that is owned by root
    beforeScript 'chmod g+w .'

	input:
    file part from ligands_parts
    file prmfile
	file protein
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

	container 'informaticsmatters/rdock-mini'
	// change permissions - see above
	beforeScript 'chmod g+w .'

	input:
	file ligands
	file part from docked_parts.collect()

	output:
	file 'results.sdf' into results

	"""
	sdsort -n -s -fSCORE docked_part*.sd |${params.score == null ? '' : " sdfilter -f'\$SCORE <= $params.score' |"}${params.nscore == null ? '' : " sdfilter -f'\$SCORE.norm <= $params.nscore' |"} sdfilter -f'\$_COUNT <= ${params.top}' > results.sdf
	"""
}

process metrics {

    beforeScript 'chmod g+w .'
    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir "$baseDir/results", mode: 'symlink'

    input:
    file 'results.sdf' from results
    file 'splitter_metrics.txt' from splitter_metrics

    output:
    file 'output.data.gz'
    file 'output.metadata'
    file 'output_metrics.txt'

    """
    python -m pipelines_utils_rdkit.filter -i results.sdf -of json -o output --meta
    mv output_metrics.txt old_metrics.txt
    grep '__InputCount__' splitter_metrics.txt | sed s/__InputCount__/DockingRDock/ > output_metrics.txt
    grep '__InputCount__' splitter_metrics.txt >> output_metrics.txt
    grep '__OutputCount__' old_metrics.txt >> output_metrics.txt
    """
}
