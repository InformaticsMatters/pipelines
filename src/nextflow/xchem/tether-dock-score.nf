#!/usr/bin/env nextflow

// tether params
params.smiles = '*.smi'
params.molfiles = '*.mol'
params.chunk_tether = 250
params.chunk_score = 1000
params.num_conformers = 1
params.ph_min = null
params.ph_max = null
params.atom_compare = 'CompareElements'
params.bond_compare = 'CompareOrder'
params.complete_rings_only = true
params.ring_matches_ring_only = true
params.minimize = 4
params.timeout_embed = null

// docking params
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.pharmafile = 'pharma.restr'
params.num_dockings = 5

// featurestein
params.fragments = 'data/mpro/hits-23.sdf.gz'

// interactions
params.iprotein = 'receptor.pdb'
params.key_hbonds = null
params.key_hydrophobic = null
params.key_halogen = null
params.key_salt_bridge = null
params.key_pi_stacking = null
params.key_pi_cation = null

// files
smilesfiles = file(params.smiles)
molfiles = file(params.molfiles)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)
pharmafile = file(params.pharmafile)
fragments = file(params.fragments)
iprotein = file(params.iprotein)


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
    publishDir '.', mode: 'copy'

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


process rdock {

    container 'informaticsmatters/rdock-mini:latest'
    errorStrategy 'retry'
    maxRetries 3

    input:
    file part from tethered_parts.flatten()
    file 'receptor.mol2' from protein
    file 'docking.prm' from prmfile
    file 'docking.as' from asfile
    file 'pharma.restr' from pharmafile

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

    publishDir ".", mode: 'copy'

	input:
    file part from docked_parts
    file fmaps

    output:
    file 'FS_*.sdf' into featurestein_parts

    """
    python -m pipelines.xchem.featurestein_score -i '$part' -if sdf -f '$fmaps' -o 'FS_${part.name[0..-4]}' -of sdf --no-gzip
    """
}

process interactions {

    container 'informaticsmatters/rdkit_pipelines:latest'
    publishDir ".", mode: 'copy'

    input:
    file part from featurestein_parts
    file iprotein

    output:
    file 'INT_*.sdf'

    """
    python -m pipelines.xchem.calc_interactions -i '$part' -if sdf -p $iprotein -o 'INT_${part.name[0..-5]}' -of sdf --no-gzip\
      ${params.key_hbond ? '--key-hbond ' + params.key_hbond : ''}\
      ${params.key_hydrophobic ? '--key-hydrophobic ' + params.key_hydrophobic : ''}\
      ${params.key_halogen ? '--key-halogen ' + params.key_halogen : ''}\
      ${params.key_salt_bridge ? '--key-salt-bridge ' + params.key_salt_bridge : ''}\
      ${params.key_pi_stacking ? '--key-pi-stacking ' + params.key_pi_stacking : ''}\
      ${params.key_pi_cation ? '--key-pi-cation ' + params.key_pi_cation : ''}
    """
}


//process xcos {
//
//    container 'informaticsmatters/rdkit_pipelines:latest'
//
//    publishDir ".", mode: 'link'
//
//	input:
//    file part from featurestein_parts
//    file fragments
//
//    output:
//    file 'XC_*.sdf.gz'
//
//    """
//    python -m pipelines.xchem.xcos -i '$part' -if sdf -f '$fragments' -o 'XC_${part.name[0..-5]}' -of sdf
//    """
//}
