#!/usr/bin/env nextflow

params.ligands = 'ligands.data.gz'
params.protein = 'protein.pdb.gz'
params.chunk = 25
params.score = null
params.limit = 0
params.digits = 4

ligands = file(params.ligands)
protein = file(params.protein)

process splitter {

    beforeScript 'chmod g+w .'

    input:
    file ligands

    output:
    file 'ligand_part*.sdf.gz' into ligand_parts mode flatten
    file 'ligand_part_metrics.txt' into splitter_metrics

    """
    python -m pipelines_utils_rdkit.filter -i $ligands -c $params.chunk -l $params.limit -d $params.digits -o ligand_part -of sdf --meta
    """
}


/* Scores each file from the ligand_parts channel sending each resulting SD file to the results channel
*/
process smog_scoring {

    beforeScript 'chmod g+w .'

	input:
    file part from ligand_parts
	file protein

    output:
    file 'scored_part*.sdf' into scored_parts

    """
	python -m pipelines.docking.smog2016 -i $part -pdb $protein -o ${part.name.replace('ligand', 'scored')[0..-8]} -of sdf --no-gzip ${params.score ? ' -t ' + params.score : ''} --threads 1 &> scored_out.log
    """
}

process joiner {

    beforeScript 'chmod g+w .'
    publishDir baseDir, pattern: "{output.data.gz,output.metadata}"

    input:
	file parts from scored_parts.collect()

	output:
	file 'output_metrics.txt' into joiner_metrics
	file 'output.data.gz'

	"""
	cat scored_part*.sdf | python -m pipelines_utils_rdkit.filter -if sdf -of json -o output --meta --thin
	"""
}

process metrics_meta {

    beforeScript 'chmod g+w .'
    publishDir baseDir


    input:
    file 'splitter_metrics.txt' from splitter_metrics
    file 'joiner_metrics.txt' from joiner_metrics

    output:
    file 'output_metrics.txt'
    file 'output.metadata'

    """
    grep '__InputCount__' splitter_metrics.txt | sed s/__InputCount__/PLI/ > output_metrics.txt
    grep '__InputCount__' splitter_metrics.txt >> output_metrics.txt
    grep '__OutputCount__' joiner_metrics.txt >> output_metrics.txt
    echo '{"type":"org.squonk.types.BasicObject","valueClassMappings":{"pliff_cscore":"java.lang.Float","pliff_iscore":"java.lang.Float","pliff_tscore":"java.lang.Float","pliff_gscore":"java.lang.Float","pliff_score":"java.lang.Float","pliff_nb_score":"java.lang.Float"}}' > output.metadata
    """
}

