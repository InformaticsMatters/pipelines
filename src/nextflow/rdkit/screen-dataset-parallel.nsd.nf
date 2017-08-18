#!/usr/bin/env nextflow

params.input = file("$baseDir/input.data.gz")
params.qsmiles
params.simmin = 0.7
params.simmax = 1.0
params.descriptor = 'rdkit'
params.metric = 'tanimoto'
params.chunk = 1000
params.limit = 0
params.digits = 4

target = file(params.input)

process splitter {

    input:
    file target

    output:
    file 'target_part*.sdf.gz' into target_parts mode flatten
    file 'target_part_metrics.txt' into splitter_metrics

    """
    python -m pipelines.rdkit.filter -i $target -c $params.chunk -l $params.limit -d $params.digits -o target_part -of sdf --meta
    """
}

process rdkitScreen {

	input:
    file part from target_parts

    output:
    file 'screened_part*.sdf.gz' into screened_parts

    """
    python -m pipelines.rdkit.screen --qsmiles '$params.qsmiles' --simmin $params.simmin --simmax $params.simmax -d $params.descriptor -m $params.metric -i $part -o ${part.name.replace('target', 'screened')[0..-8]} -of sdf
    """
}

process joiner {

    publishDir baseDir, pattern: "{output.data.gz,output.metadata}"

    input:
	file parts from screened_parts.collect()

	output:
	file 'output_metrics.txt' into joiner_metrics
	file 'output.data.gz'
	file 'output.metadata'

	"""
	zcat $parts | python -m pipelines.rdkit.filter -if sdf -of json -o output --meta
	"""
}

process metrics {

    publishDir baseDir

    input:
    file 'splitter_metrics.txt' from splitter_metrics
    file 'joiner_metrics.txt' from joiner_metrics

    output:
    file 'output_metrics.txt'

    """
    grep '__InputCount__' splitter_metrics.txt | sed s/__InputCount__/RDKitScreen/ > output_metrics.txt
    grep '__InputCount__' splitter_metrics.txt >> output_metrics.txt
    grep '__OutputCount__' joiner_metrics.txt >> output_metrics.txt
    """
}

