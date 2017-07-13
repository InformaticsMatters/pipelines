#!/usr/bin/env nextflow

/* Nextflow pipline that runs MolAlign from Ling Chan <slchan.bc@gmail.com>
* This is designed to be run in the informaticsmatters/rdkit_pipelines_molalign container.
*/


params.conformers = 'input.data.gz'
params.field = 'StructureNum'
params.assemblies = 5


conformers = file(params.conformers)

/* Splits the input SD file into multiple files corresponding to each input molecule.
* Assumes that the input contains one or more conformers of multiple input structures and
* a field named ${params.field} identifies the ID of the source structure so that the
* conformers can be split into separate files for each structure.
*/
process split {

	input:
    file conformers

    output:
    file 'molecule_*.sdf' into conformer_parts
    file 'molecule__metrics.txt' into input_metrics


    """
    python -m pipelines.rdkit.splitter -i $conformers -f $params.field -o molecule_ --meta
    """
}

process molalign {

	input:
    file parts from conformer_parts

    output:
    file 'A.sdf' into assemblies

    """
    echo "$parts" | sed "s/ /\\n/g" > InpFile.txt
    echo "maxout $params.assemblies" > params.txt
    MolAlignA params=params.txt outdir=. inpconfs=InpFile.txt
    """
}

process results {

	publishDir './', mode: 'copy'

	input:
    file assemblies
    file 'input_metrics.txt' from input_metrics

    output:
    file 'output.data.gz'
    file 'output.metadata'
    file 'output_metrics.txt'

    """
    python -m pipelines.rdkit.filter -i $assemblies -if sdf -o output -of json --rename uuid:ConformerUUID
    echo '{"type": "org.squonk.types.MoleculeObject"}' > output.metadata
    grep '__InputCount__' input_metrics.txt > output_metrics.txt
    count=`fgrep -c '\$\$\$\$' $assemblies`
    echo "__OutputCount__=\$count" >> output_metrics.txt
    echo "MolAlignAssemblies=\$count" >> output_metrics.txt
    """
}


