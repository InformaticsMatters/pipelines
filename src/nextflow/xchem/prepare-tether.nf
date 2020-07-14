#!/usr/bin/env nextflow

params.candidates = "data/mpro/expanded-17.json"
params.chunk = 1000
params.limit = 0
params.digits = 4
params.generate_filenames = false

candidates = file(params.candidates)

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
    file smiles from smiles.flatten() //collect().toSortedList().flatten()
    file mol from mols.flatten() //collect().toSortedList().flatten()

    output:
    file 'Tethered_*.sdf'
    """
    python -m pipelines.xchem.prepare_tether --smi '$smiles' --mol '$mol' -o 'Tethered_${smiles.name[0..-5]}'
    """
}
