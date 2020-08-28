#!/usr/bin/env nextflow

params.smiles = '*.smi'
params.molfiles = '*.mol'
params.fragments = "data/mpro/hits-17.sdf.gz"
params.chunk_tether = 250
params.chunk_score = 10000
params.limit = 0
params.digits = 4
params.generate_filenames = false
params.num_conformers = 10

// files
smilesfiles = file(params.smiles)
molfiles = file(params.molfiles)
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
