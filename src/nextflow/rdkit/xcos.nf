#!/usr/bin/env nextflow

params.inputs = "data/mpro/poses.sdf"
params.fragments = "data/mpro/hits-17.sdf"
params.chunk = 500
params.limit = 0
params.digits = 4

inputs = file(params.inputs)
fragments = file(params.fragments)

process splitter {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file inputs

    output:
    file 'inputs_part*.sdf.gz' into inputs_parts mode flatten

    """
    python -m pipelines_utils_rdkit.filter -i '$inputs' -c $params.chunk -l $params.limit -d $params.digits -o 'inputs_part_' -of sdf
    """
}

process xcos {

    container 'informaticsmatters/rdkit_pipelines:latest'

	input:
    file part from inputs_parts
    file fragments

    output:
    file 'scored_part*.sdf' into scored_parts

    """
    python -m pipelines.rdkit.xcos -i '$part' -f '$fragments' -o '${part.name.replace('inputs', 'scored')[0..-8]}' -of sdf --no-gzip
    """
}

process joiner {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir ".", mode: 'move'

    input:
	file parts from scored_parts.collect()

	output:
	file 'xcos_scored.sdf.gz'

	"""
	cat '$parts' | gzip > xcos_scored.sdf.gz
	"""
}
