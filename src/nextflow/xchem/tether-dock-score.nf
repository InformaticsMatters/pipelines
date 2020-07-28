#!/usr/bin/env nextflow

// splitter params
params.chunk_size_expand = 200
params.limit = 0
params.digits = 4

// tether params
params.smiles = '*.smi'
params.molfiles = '*.mol'
params.ph_min = null
params.ph_max = null
params.timeout_embed = null
params.chunk_size_tether =  500

// docking params
params.protein = 'data/mpro/Mpro-x0387_0.mol2'
params.prmfile = 'data/mpro/docking-tethered.prm'
params.asfile = 'data/mpro/docking-tethered.as'
params.num_dockings = 5

// featurestein
params.fragments = 'data/mpro/hits-23.sdf.gz'

// files
smiles = file(params.smiles)
molfiles = file(params.molfiles)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)
fragments = file(params.fragments)


process splitter {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file smiles from smiles.flatten()
    file mol from molfiles.flatten()

    output:
    file '*.mol' into split_mols
    file '*.smi' into split_smiles

    """
    stem=${smiles.name[0..-5]}
    split -l $params.chunk_size_expand -d -a 3 --additional-suffix .smi $smiles \${stem}_
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
    file mol from split_mols.flatten()
    file smiles from split_smiles.flatten()

    output:
    file 'Tethered_*.sdf' into tethered_parts

    """
    python -m pipelines.xchem.prepare_tether --smi '$smiles' --mol '$mol'\
      --chunk-size $params.chunk_size_tether\
      ${params.ph_min != null ? '--min-ph ' + params.ph_min : ''}\
      ${params.ph_max != null ? '--max-ph ' + params.ph_max : ''}\
      ${params.timeout_embed != null ? '--timeout-embed ' + params.timeout_embed : ''}\
      -o 'Tethered_${smiles.name[0..-5]}'
    """
}

process rdock {

    container 'informaticsmatters/rdock-mini:latest'
    errorStrategy 'retry'
    maxRetries 3

	input:
    file part from tethered_parts.flatten()
	file 'receptor.mol2' from protein
	file 'docking.prm' from prmfile
	file 'docking.as' from asfile

    output:
    file 'Docked_*.sd' optional true into docked_parts

    """
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('Tethered', 'Docked')[0..-5]} > docked_out.log
    """
}

process gen_feat_maps {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file fragments

    output:
    file 'featurestein.p' into fmaps

    """
    python -m pipelines.xchem.featurestein_generate -i '$fragments' -f featurestein.p
    """
}

process featurestein {

    container 'informaticsmatters/rdkit_pipelines:latest'

	input:
    file part from docked_parts
    file fmaps

    output:
    file 'FS_*.sdf' into featurestein_parts

    """
    python -m pipelines.xchem.featurestein_score -i '$part' -if sdf -f '$fmaps' -o 'FS_${part.name[0..-4]}' -of sdf --no-gzip
    """
}

process xcos {

    container 'informaticsmatters/rdkit_pipelines:latest'

    publishDir ".", mode: 'link'

	input:
    file part from featurestein_parts
    file fragments

    output:
    file 'XC_*.sdf.gz'

    """
    python -m pipelines.xchem.xcos -i '$part' -if sdf -f '$fragments' -o 'XC_${part.name[0..-5]}' -of sdf
    """
}
