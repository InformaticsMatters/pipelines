#!/usr/bin/env nextflow

params.inputs = "data/mpro/poses.sdf.gz"
params.fragments = "data/mpro/hits-17.sdf.gz"
params.chunk = 5000
params.limit = 0
params.digits = 4

inputs = file(params.inputs)
fragments = file(params.fragments)

process generate_feat_maps {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file fragments

    output:
    file 'featurestein.p' into fmaps

    """
    python -m pipelines.xchem.featurestein_generate -i '$fragments' -f featurestein.p
    """
}

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

process score {

    container 'informaticsmatters/rdkit_pipelines:latest'

	input:
    file part from inputs_parts
    file fmaps

    output:
    file 'scored_part*.sdf' into scored_parts

    """
    python -m pipelines.xchem.featurestein_score -i '$part' -f '$fmaps' -o '${part.name.replace('inputs', 'scored')[0..-8]}' -of sdf --no-gzip
    """
}

process joiner {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir ".", mode: 'link'

    input:
	file parts from scored_parts.collect()

	output:
	file 'featurestein_scored.sdf.gz'

	"""
	cat $parts | gzip > featurestein_scored.sdf.gz
	"""
}
