#!/usr/bin/env nextflow

/* Squonk Nextflow pipline that runs Docking using rDock, filtering the poses relative to the score from docking the
* reference ligand.
*
* To test this manually run something like this:
* nextflow -c src/nextflow/nextflow-docker.config run src/nextflow/docking/rdock-filter.nsd.nf --ligands data/nudt7/ligands.data.gz --refmol data/nudt7/refmol.mol --receptor data/nudt7/receptor.mol2 --num_dockings 5
*/

params.refmol = "$baseDir/refmol.mol"
params.ligands = "$baseDir/ligands.data.gz"
params.receptor = "$baseDir/receptor.mol2.gz"
params.chunk = 25
params.num_dockings = 25
params.top = 1
params.limit = 0
params.digits = 4
params.threshold = 0.0
params.field = 'SCORE.norm'

refmol = file(params.refmol)
ligands = file(params.ligands)
receptor = file(params.receptor)

expl2 = Channel.value( "Hello there $receptor" )

process create_cavity {

    container 'informaticsmatters/rdock-mini:latest'
    beforeScript 'chmod g+w .'

    input:
    file refmol
    file receptor

    output:
    file 'receptor.prm' into prmfile
    file 'receptor.as' into asfile

    """
    gunzip -c $receptor > receptor.mol2
    cat << EOF > receptor.prm
RBT_PARAMETER_FILE_V1.00
RECEPTOR_FILE receptor.mol2
RECEPTOR_FLEX 3.0
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL $refmol
    RADIUS 3.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
EOF

rbcavity -was -d -r receptor.prm > rbcavity.log
    """
}

/* Docks the reference ligand
*/
process dock_reference_ligand {

    container 'informaticsmatters/rdock-mini:latest'
    beforeScript 'chmod g+w .'

    publishDir "$baseDir/results", mode: 'copy'

    input:
    file receptor
    file 'receptor.as' from asfile
    file 'receptor.prm' from prmfile
    file 'refmol.mol' from refmol

    output:
    file 'best_ligand.sdf' into best_ligand

    """
    gunzip -c $receptor > receptor.mol2
    rbdock -i refmol.mol -r receptor.prm -p dock.prm -n $params.num_dockings -o docked_ligand > docked_ligand_out.log
    sdsort -n -s -fSCORE docked_ligand.sd | sdfilter -f'\$_COUNT <= 1' > best_ligand.sdf
    """
}

/* Splits the input into multiple files of ${params.chunk} records.
*/
process splitter {

    //beforeScript 'chmod g+w .'
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

/* Docks each file from the ligands_parts channel sending each resulting SD file to the results channel
*/
process dock_ligands {

    container 'informaticsmatters/rdock-mini:latest'
    // change permissions on the work dir so that the rdock user in the container
    // can write to the directory that is owned by root
    beforeScript 'chmod g+w .'

    input:
    file part from ligands_parts
    file receptor
    file 'receptor.as' from asfile
    file 'receptor.prm' from prmfile

    output:
    file 'docked_part*.sd' into docked_parts

    """
    gunzip -c $receptor > receptor.mol2
    rbdock -i $part -r receptor.prm -p dock.prm -n $params.num_dockings -o ${part.name.replace('ligands', 'docked')[0..-5]} > docked_out.log
    """
}

/* Filter, combine and publish the results.
* Poses are only included if they are within ${params.threshold} of the best score obtained from docking the
* reference ligand into the same receptor (output of the dock_ligand process).
*/
process combine_and_filter {

	container 'informaticsmatters/rdock-mini:latest'
	beforeScript 'chmod g+w .'

	input:
	file parts from docked_parts.collect()
	file best from best_ligand

	output:
	file 'rdock_results.sdf' into results

	"""
	FSCORE=\$(sdreport -nh -t${params.field} best_ligand.sdf | cut -f 2 | awk '{\$1=\$1};1')
	ASCORE=\$(awk "BEGIN {print \$FSCORE + ${params.threshold}}")
	echo "Processing $parts with normalised score filter of \$ASCORE"
	sdsort -n -s -f${params.field} docked_part*.sd | sdfilter -f"\\\$${params.field} < \$ASCORE" | sdfilter -f'\$_COUNT <= ${params.top}' > rdock_results.sdf
	"""
}

process results {

    beforeScript 'chmod g+w .'
    container 'informaticsmatters/rdkit_pipelines:latest'
    beforeScript 'chmod g+w .'

    publishDir "$baseDir/results", mode: 'copy'

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
    echo -n 'DockingRDock=' >> output_metrics.txt
    echo \$((`grep '__InputCount__' splitter_metrics.txt | cut -d '=' -f 2` * ${params.num_dockings})) >> output_metrics.txt
    grep '__InputCount__' splitter_metrics.txt >> output_metrics.txt
    grep '__OutputCount__' old_metrics.txt >> output_metrics.txt
    """
}
