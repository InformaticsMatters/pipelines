#!/usr/bin/env nextflow

/* Nextflow pipline that runs MolAlign from Ling Chan <slchan.bc@gmail.com>
*/


params.conformers = 'input.data.gz'
params.field = 'UUID'
params.max_overlays = 5


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
    
    
    """
    python -m pipelines.rdkit.splitter -i $conformers -f $params.field -o molecule_
    """
}

process molalign {

	input:
    file parts from conformer_parts
   
    """
    echo "$parts" | sed "s/ /\\n/g" > InpFile.txt
    echo "maxout $params.max_overlays" > params.txt
    MolAlign params=params.txt outdir=results inpconfs=InpFile.txt
    cat results/A*.sdf > all.sdf
    """

}
