from __future__ import print_function
import sys, gzip
from rdkit import Chem
from sanifix import fix_mol
from ijson import items


'''
Log output to STDERR
'''
def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def add_default_io_args(parser):
    parser.add_argument('-i', '--input', help="Input SD file, if not defined the STDIN is used")
    parser.add_argument('-o', '--output', help="Base name for output file (no extension). If not defined then SDTOUT is used for the structures and output is used as base name of the other files.")
    parser.add_argument('-if', '--informat', choices=['sdf', 'json'], help="Input format. When using STDIN this must be specified.")
    #parser.add_argument('-of', '--outformat', choices=['sdf', 'json'], help="Output format. When using STDOUT this must be specified.")


'''
Default approach to handling the inputs and outputs
'''
def default_open_input_output(inputDef, inputFormat, outputDef, defaultOutput):
    input, suppl = default_open_input(inputDef, inputFormat)

    if outputDef:
        output = gzip.open(outputDef + '.sdf.gz','w+')
        outputBase = outputDef
    else:
        output = sys.stdout
        outputBase = defaultOutput
        
    writer = Chem.SDWriter(output)
    
    return input,output,suppl,writer,outputBase


def default_open_input(inputDef, inputFormat):
    if not inputDef and not inputFormat:
        raise ValueError('Must specify either an input file name or an input format (or both)')
    elif inputFormat == 'sdf' or (inputDef and (inputDef.lower().endswith('.sdf') or inputDef.lower().endswith('.sdf.gz'))):
        input, suppl = default_open_sdf(inputDef)
    elif inputFormat == 'json' or (inputDef and (inputDef.lower().endswith('.data') or inputDef.lower().endswith('.data.gz'))):
        input, suppl = default_open_json(inputDef)
    else:
        raise ValueError('Unsupported input format')

    return input, suppl


def default_open_sdf(inputDef):
    if inputDef:
        if inputDef.lower().endswith('.gz'):
            input = gzip.open(inputDef)
        else:
            input = open(inputDef, 'r')
    else:
        input = sys.stdin
    suppl = Chem.ForwardSDMolSupplier(input)
    return input, suppl


def default_open_json(inputDef):
    if inputDef:
        if inputDef.lower().endswith('.gz'):
            input = gzip.open(inputDef)
        else:
            input = open(inputDef, 'r')
    else:
        input = sys.stdin
    suppl = generate_mols_from_json(input)
    return input, suppl


'''
Write the metrics data
'''
def write_metrics(baseName, values):
    m = open(baseName  + '_metrics.txt', 'w')
    for key in values:
        m.write(key + '=' + str(values[key]) + "\n")
    m.flush()
    m.close()


def parse_mol_simple(my_type, txt):
    """Function to parse individual mols given a type"""
    if my_type == "mol":
        # Try this way
        mol = Chem.MolFromMolBlock(txt.strip())
        if mol is None:
            mol = Chem.MolFromMolBlock(txt)
        if mol is None:
            mol = Chem.MolFromMolBlock("\n".join(txt.split("\n")[1:]))
        # Now try to do sanidfix
        if mol is None:
            mol = fix_mol(Chem.MolFromMolBlock(txt, False))
        # And again
        if mol is None:
            mol = fix_mol(Chem.MolFromMolBlock(txt.strip(), False))
    elif my_type == "smiles":
        # Assumes that smiles is the first column -> and splits on chemaxon
        mol = Chem.MolFromSmiles(txt.split()[0].split(":")[0])
    if mol is None:
        log('Failed to parse mol', txt)
    return mol


'''
Function to get the RDKit mol from MoleculeObject JSON
'''
def parse_mol_json(molobj):
    #print "reading mol",molobj["uuid"],molobj["format"]
    molstr = str(molobj["source"])
    # Get the format and use this as a starting point to work out
    molformat = molobj["format"]
    # Now parse it with RDKit
    mol = parse_mol_simple(molformat, molstr)
    mol.SetProp("uuid", str(molobj["uuid"]))
    values = molobj["values"]
    if values:
        for key in values:
            mol.SetProp(str(key), str(values[key]))
    return mol


'''
Create a supplier of RDKit Mol objects from the json
'''
def generate_mols_from_json(json):
    j=0
    for item in items(json, "item"):
        j+=1
        #print "  reading",j
        mol = parse_mol_json(item)
        yield mol