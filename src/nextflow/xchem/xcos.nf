#!/usr/bin/env nextflow

params.inputs = 'data/mpro/poses.sdf.gz'
params.fragments = 'data/mpro/hits-17.sdf.gz'
params.threshold = 0.4 // XCos score theshold
params.chunk = 500     // chunk size to split input into
params.limit = 0       // max number of molecules to process
params.digits = 4      // number of digits for the split file name number

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
    python -m pipelines.xchem.xcos -i '$part' -f '$fragments' -t $params.threshold -o '${part.name.replace('inputs', 'scored')[0..-8]}' -of sdf --no-gzip
    """
}

process joiner {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir ".", mode: 'link'

    input:
	file parts from scored_parts.collect()

	output:
	file 'xcos_scored.sdf.gz'

	"""
	cat $parts | gzip > xcos_scored.sdf.gz
	"""
}
