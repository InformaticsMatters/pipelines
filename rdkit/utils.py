from __future__ import print_function
import sys, gzip
from rdkit import Chem

### start function defintions #########################################
          
def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    
def defaultOpenInputOutput(inputDef, outputDef, defaultOutput):
	if inputDef:
		input = gzip.open(inputDef)
	else:
		input = sys.stdin
	suppl = Chem.ForwardSDMolSupplier(input)

	if outputDef:
		output = gzip.open(outputDef + '.sdf.gz','w+')
		output_base = outputDef
	else:
		output = sys.stdout
		output_base = defaultOutput
		
	writer = Chem.SDWriter(output)
	
	return input,output,suppl,writer,output_base
