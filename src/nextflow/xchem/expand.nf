#!/usr/bin/env nextflow

// expand params
params.hits = 'data/mpro/hits-5.sdf.gz'
params.token = null
params.hac_min = 3
params.hac_max = 3
params.rac_min = 1
params.rac_max = 1
params.hops = 1
params.server = null


// files
hits = file(params.hits)

process fragnet_expand {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir ".", mode: 'copy'

    input:
    file hits

    output:
    file '*.smi' into smiles
    file '*.mol' into mols

    """
    python -m pipelines.xchem.fragnet_expand -i '$hits' ${params.token ? '--token ' + params.token : ''}\
      --hops $params.hops\
      --hac-min $params.hac_min\
      --hac-max $params.hac_max\
      --rac-min $params.rac_min\
      --rac-max $params.rac_max\
      ${params.server ? '--server ' + params.server : ''}
    """
}