#!/usr/bin/env nextflow

params.candidates = "data/mpro/expanded-17.json"
params.fragments = "data/mpro/hits-17.sdf.gz"
params.chunk_tether = 250
params.chunk_score = 10000
params.limit = 0
params.digits = 4
params.generate_filenames = false
params.num_conformers = 10

candidates = file(params.candidates)
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

process split_json {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file candidates

    output:
    file '*.smi' into smiles
    file '*.mol' into mols

    """
    python -m pipelines.xchem.split_fragnet_candidates -i '$candidates' ${params.generate_filenames ? '--generate-filenames' : ''}
    """
}

process tether {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file smiles from smiles.flatten()
    file mol from mols.flatten()

    output:
    file 'Tethered_*.sdf' into tethered_parts

    """
    python -m pipelines.xchem.prepare_tether --smi '$smiles' --mol '$mol' --chunk-size $params.chunk_score --num-conformers $params.num_conformers -o 'Tethered_${smiles.name[0..-5]}'
    """
}

process score {

    container 'informaticsmatters/rdkit_pipelines:latest'
    publishDir '.'

	input:
    file part from tethered_parts.flatten()
    file fmaps

    output:
    file 'Scored_*.sdf' into scored_parts

    """
    python -m pipelines.xchem.featurestein_score -i '$part' -f '$fmaps' -o 'Scored_${part.name[0..-5]}' -of sdf --no-gzip
    """
}
