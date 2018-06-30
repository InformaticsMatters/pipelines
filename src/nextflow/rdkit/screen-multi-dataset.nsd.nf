#!/usr/bin/env nextflow

params.target = "$baseDir/target.data.gz"
params.query = "$baseDir/query.data.gz"
params.simmin = 0.7
params.simmax = 1.0
params.descriptor = 'rdkit'
params.metric = 'tanimoto'
params.chunk = 500
params.limit = 0
params.digits = 4

target = file(params.target)
query = file(params.query)

process splitter {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file target

    output:
    file 'target_part*.sdf.gz' into target_parts mode flatten
    file 'target_part_metrics.txt' into splitter_metrics

    """
    python -m pipelines_utils_rdkit.filter -i $target -c $params.chunk -l $params.limit -d $params.digits -o target_part -of sdf --meta
    """
}

process rdkitScreen {

    container 'informaticsmatters/rdkit_pipelines:latest'

	input:
    file part from target_parts
    file 'query.data.gz' from query

    output:
    file 'screened_part*.sdf.gz' into screened_parts

    """
    python -m pipelines.rdkit.screen_multi -i $part --qjson query.data.gz --simmin $params.simmin --simmax $params.simmax -d $params.descriptor -m $params.metric -o ${part.name.replace('target', 'screened')[0..-8]} -of sdf
    """
}

process joiner {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir baseDir, mode: 'link'

    input:
    file 'splitter_metrics.txt' from splitter_metrics
	file parts from screened_parts.collect()

	output:
	file 'output_metrics.txt'
	file 'output.data.gz'
	file 'output.metadata'

	"""
	zcat $parts | python -m pipelines_utils_rdkit.filter -if sdf -of json -o output --meta
	mv output_metrics.txt joiner_metrics.txt
	grep '__InputCount__' splitter_metrics.txt | sed s/__InputCount__/RDKitScreen/ > output_metrics.txt
    grep '__InputCount__' splitter_metrics.txt >> output_metrics.txt
    grep '__OutputCount__' joiner_metrics.txt >> output_metrics.txt
	"""
}
