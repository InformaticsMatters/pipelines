#!/usr/bin/env nextflow

/* Nextflow pipline that runs MolAlign1o1 from Ling Chan <slchan.bc@gmail.com>
*/


params.conformers = 'input.data.gz'
params.template = 'first_mol.sdf'
params.field = 'StructureNum'
params.nout = null


conformers = file(params.conformers)

template = file(params.template)

/* Splits the input SD file into multiple files corresponding to each input molecule.
* Assumes that the input contains one or more conformers of multiple input structures and 
* a field named ${params.field} identifies the ID of the source structure so that the
* conformers can be split into separate files for each structure. 
*/
process split {

	container 'informaticsmatters/rdkit_pipelines'

	input:
    file conformers

    output:
    file 'molecule_*.sdf' into conformer_parts
    
    
    """
    python -m pipelines.rdkit.splitter -i $conformers -f $params.field -o molecule_
    """
}

/* Performs the alignment for the conformers of a single input structure
*
*/
process molalign {

	container 'informaticsmatters/molalign'

	input:
	file template
    file parts from conformer_parts
    
    output:
    file 'result_*.sdf' into assemblies
   
   	"""
   	for f in molecule_*.sdf
   	do
   	    name=\${f/molecule_/result_}
   	    MolAlign1o1 $template \$f ${params.nout ? 'nout=' + params.nout : ''} outsdf=\$name outscore=\${name%.sdf}.sco
   	done
   	"""
}

/* Collects the resulting SD files and converts back to JSON
* TODO - generate metadata
*
*/
process results {

	container 'informaticsmatters/rdkit_pipelines'

	publishDir './', mode: 'copy'
	
	input:
	file results from assemblies
	
	output:
	file 'results.data.gz'
	

	"""
	cat $results | python -m pipelines.rdkit.filter -if sdf -o results -of json
	"""
}

