#!/usr/bin/env nextflow

/* Example Nextflow pipline that runs SMoG2016 scoring
*/


params.ligands = 'ligands.sdf.gz'
params.protein = 'protein.pdb'
params.chunk = 25
params.score = null


ligands = file(params.ligands)
protein = file(params.protein)

/* Splits the input SD file into multiple files of ${params.chunk} records.
* Each file is sent individually to the ligand_parts channel.
* Parts are renamed so as to be in correct sorted area.
*/
process sdsplit {

	input:
    file ligands

    output:
    file 'ligands_part*' into ligand_parts mode flatten
    
    
    """
    python -m pipelines_utils_rdkit.filter -i $ligands -c $params.chunk -d 5 -o ligands_part -of sdf
    """
}

/* Scores each file from the ligand_parts channel sending each resulting SD file to the results channel
*/
process smog_scoring {

	input:
    file part from ligand_parts
	file protein

    output:
    file 'scored_part*.sdf' into scored_parts
    
    """
	python -m pipelines.docking.smog2016 -i $part -pdb $protein -o ${part.name.replace('ligands', 'scored')[0..-8]} -of sdf --thin --no-gzip ${params.score ? ' -t ' + params.score : ''} --threads 1 &> scored_out.log
    """
 }

/* Recombine the results
*/
process results {

	
	input:
	file ligands
	file part from scored_parts.collect()
	
	output:
	file 'output.sdf.gz' into results
	
 
	"""
	cat scored_part*.sdf | gzip > output.sdf.gz
	"""
}

results.println { "Results: $it" }