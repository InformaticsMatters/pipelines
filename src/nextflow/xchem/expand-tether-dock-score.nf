#!/usr/bin/env nextflow

// expand params
params.hits = 'data/mpro/hits-17.sdf.gz'
params.chunk_size_expand = 200
params.limit = 0
params.digits = 4
params.token = null
params.hac_min = 3
params.hac_max = 3
params.rac_min = 1
params.rac_max = 1
params.hops = 1
params.server = null

// tether params
params.ph_min = null
params.ph_max = null
params.timeout_embed = null
params.chunk_size_tether =  null

// docking params
params.protein = 'data/mpro/Mpro-x0387_0.mol2'
params.prmfile = 'data/mpro/docking-tethered.prm'
params.asfile = 'data/mpro/docking-tethered.as'
params.num_dockings = 5

// featurestein

// files
hits = file(params.hits)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)

process fragnet_expand {

    container 'informaticsmatters/rdkit_pipelines:latest'

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

process splitter {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file smiles from smiles.flatten()
    file mol from mols.flatten()

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
    file 'Docked_*.sd' into docked_parts

    """
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('Tethered', 'Docked')[0..-5]} > docked_out.log
    """
}

process gen_feat_maps {

    container 'informaticsmatters/rdkit_pipelines:latest'

    input:
    file hits

    output:
    file 'featurestein.p' into fmaps

    """
    python -m pipelines.xchem.featurestein_generate -i '$hits' -f featurestein.p
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

    publishDir ".", mode: 'move'

	input:
    file part from featurestein_parts
    file hits

    output:
    file 'XC_*.sdf'

    """
    python -m pipelines.xchem.xcos -i '$part' -if sdf -f '$hits' -o 'XC_${part.name[0..-5]}' -of sdf --no-gzip
    """
}
