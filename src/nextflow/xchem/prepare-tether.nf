#!/usr/bin/env nextflow

params.smiles = '*.smi'
params.molfiles = '*.mol'
params.chunk_tether = 250
params.chunk_score = 1000
params.limit = 0
params.num_conformers = 1
params.ph_min = null
params.ph_max = null
params.atom_compare = 'CompareElements'
params.bond_compare = 'CompareOrder'
params.complete_rings_only = true
params.ring_matches_ring_only = true
params.minimize = 4
params.timeout_embed = null


smilesfiles = file(params.smiles)
molfiles = file(params.molfiles)

process splitter {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file smiles from smilesfiles.flatten()
    file mol from molfiles.flatten()

    output:
    file '*.mol' into mols
    file '*.smi' into smiles

    """
    stem=${smiles.name[0..-5]}
    split -l $params.chunk_tether -d -a 3 --additional-suffix .smi $smiles \${stem}_
    mv $smiles ${smiles}.orig
    for f in *.smi
    do
      cp $mol \${f:0:-4}.mol
    done
    mv $mol ${mol}.orig
    """
}

process tether {

    container 'informaticsmatters/rdkit_pipelines:latest'
    publishDir '.'

    input:
    file smiles from smiles.flatten()
    file mol from mols.flatten()

    output:
    file 'Tethered_*.sdf' into tethered_parts

    """
    python -m pipelines.xchem.prepare_tether --smi '$smiles' --mol '$mol' --chunk-size $params.chunk_score\
      -o 'Tethered_${smiles.name[0..-5]}'\
      --num-conformers $params.num_conformers\
      --atom-compare $params.atom_compare --bond-compare $params.bond_compare\
      ${params.complete_rings_only ? '--complete-rings-only' : ''}\
      ${params.ring_matches_ring_only ? '--ring-matches-ring-only' : ''}\
      --minimize $params.minimize\
      ${params.ph_min != null ? '--min-ph ' + params.ph_min : ''}\
      ${params.ph_max != null ? '--max-ph ' + params.ph_max : ''}\
      ${params.timeout_embed != null ? '--timeout-embed ' + params.timeout_embed : ''}
    """
}
