#!/usr/bin/env nextflow

params.chunk = 25
params.limit = 0
params.digits = 4

// docking params
params.ligands = 'ligands.sdf'
params.protein = 'receptor.mol2'
params.prmfile = 'docking.prm'
params.asfile = 'docking.as'
params.pharmafile = 'pharma.restr'
params.num_dockings = 5

// featurestein
params.fragments = 'data/mpro/hits-23.sdf.gz'

// interactions
params.iprotein = 'receptor.pdb'
params.key_hbond = null
params.key_hydrophobic = null
params.key_halogen = null
params.key_salt_bridge = null
params.key_pi_stacking = null
params.key_pi_cation = null

// files
ligands = file(params.ligands)
protein = file(params.protein)
prmfile = file(params.prmfile)
asfile = file(params.asfile)
pharmafile = file(params.pharmafile)
fragments = file(params.fragments)
iprotein = file(params.iprotein)

process sdsplit {

    container 'informaticsmatters/rdock-mini:latest'

    input:
    file ligands

    output:
    file 'ligands_part*.sd' into ligand_parts

    """
    sdsplit -${params.chunk} -oligands_part_ $ligands

    if [ -f ligands_part_1.sd ]; then mv ligands_part_1.sd ligands_part_01.sd; fi
    if [ -f ligands_part_2.sd ]; then mv ligands_part_2.sd ligands_part_02.sd; fi
    if [ -f ligands_part_3.sd ]; then mv ligands_part_3.sd ligands_part_03.sd; fi
    if [ -f ligands_part_4.sd ]; then mv ligands_part_4.sd ligands_part_04.sd; fi
    if [ -f ligands_part_5.sd ]; then mv ligands_part_5.sd ligands_part_05.sd; fi
    if [ -f ligands_part_6.sd ]; then mv ligands_part_6.sd ligands_part_06.sd; fi
    if [ -f ligands_part_7.sd ]; then mv ligands_part_7.sd ligands_part_07.sd; fi
    if [ -f ligands_part_8.sd ]; then mv ligands_part_8.sd ligands_part_08.sd; fi
    if [ -f ligands_part_9.sd ]; then mv ligands_part_9.sd ligands_part_09.sd; fi
    """
}

process rdock {

    container 'informaticsmatters/rdock-mini:latest'
    errorStrategy 'retry'
    maxRetries 3

    input:
    file part from ligand_parts.flatten()
    file 'receptor.mol2' from protein
    file 'docking.prm' from prmfile
    file 'docking.as' from asfile
    file 'pharma.restr' from pharmafile

    output:
    file 'Docked_*.sd' optional true into docked_parts

    """
    rbdock -r $prmfile -p dock.prm -n $params.num_dockings -i $part -o ${part.name.replace('ligands', 'Docked')[0..-4]} > docked_out.log
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
