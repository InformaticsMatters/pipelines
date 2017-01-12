from __future__ import print_function
import sys, gzip, json, uuid
from rdkit import Chem
from sanifix import fix_mol
from StreamJsonListLoader import StreamJsonListLoader



def log(*args, **kwargs):
    """
    Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def add_default_io_args(parser):
    parser.add_argument('-i', '--input', help="Input SD file, if not defined the STDIN is used")
    parser.add_argument('-o', '--output', help="Base name for output file (no extension). If not defined then SDTOUT is used for the structures and output is used as base name of the other files.")
    parser.add_argument('-if', '--informat', choices=['sdf', 'json'], help="Input format. When using STDIN this must be specified.")
    parser.add_argument('-of', '--outformat', choices=['sdf', 'json'], help="Output format. Defaults to 'sdf'.")
    parser.add_argument('--meta', action='store_true', help='Write metadata and metrics files')


def default_open_input_output(inputDef, inputFormat, outputDef, defaultOutput, outputFormat):
    """Default approach to handling the inputs and outputs"""
    input, suppl = default_open_input(inputDef, inputFormat)
    output,writer,outputBase = default_open_output(outputDef, defaultOutput, outputFormat)
    return input,output,suppl,writer,outputBase


def default_open_input(inputDef, inputFormat):
    if not inputDef and not inputFormat:
        raise ValueError('Must specify either an input file name or an input format (or both)')
    elif inputFormat == 'sdf' or (inputDef and (inputDef.lower().endswith('.sdf') or inputDef.lower().endswith('.sdf.gz'))):
        input, suppl = default_open_input_sdf(inputDef)
    elif inputFormat == 'json' or (inputDef and (inputDef.lower().endswith('.data') or inputDef.lower().endswith('.data.gz'))):
        input, suppl = default_open_input_json(inputDef)
    else:
        raise ValueError('Unsupported input format')

    return input, suppl


def default_open_input_sdf(inputDef):
    """Open the input as a SD file (possibly gzipped if ending with .gz) according to RDKit's ForwardSDMolSupplier

    :param inputDef: The name of the file. If None then STDIN is used (and assumed not to be gzipped)
    """
    if inputDef:
        input = open_file(inputDef)
    else:
        input = sys.stdin
    suppl = Chem.ForwardSDMolSupplier(input)
    return input, suppl


def default_open_input_smiles(inputDef, delimiter='\t', smilesColumn=0, nameColumn=1, titleLine=False):
    """Open the input as a file of smiles (possibly gzipped if ending with .gz) according to RDKit's SmilesMolSupplier

    :param inputDef: The name of the file. If None then STDIN is used (and assumed not to be gzipped)
    """
    if inputDef:
        input = open_file(inputDef)
    else:
        input = sys.stdin
    # SmilesMolSupplier is a bit strange as it can't accept a file like object!
    txt = input.read()
    input.close()
    suppl = Chem.SmilesMolSupplier()
    suppl.SetData(txt, delimiter=delimiter, smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine)
    return suppl


def open_file(filename):
    """Open the file as a SDF gunzipping it if it ends with .gz"""
    if filename.lower().endswith('.gz'):
        return gzip.open(filename)
    else:
        return open(filename, 'r')

def open_smarts(filename):
    """Very simple smarts parser that expects smarts expression (and nothing else) on each line (no header)"""
    f = open(filename)
    count = 0
    for line in f:
        count += 1
        mol = Chem.MolFromSmarts(line)
        if mol:
            mol.SetIntProp("idx", count)
        yield mol


def default_open_input_json(inputDef, lazy=True):
    """Open the given input as JSON array of Squonk MoleculeObjects

    :param inputDef: The name of the input file, or None if to use STDIN. If filename ends with .gz will be gunzipped
    :param lazy: Use lazy loading of the JSON. If True will allow handling of large datasets without being loaded into memory,
    but may be less robust and will be slower.
    """
    if inputDef:
        if inputDef.lower().endswith('.gz'):
            input = gzip.open(inputDef)
        else:
            input = open(inputDef, 'r')
    else:
        input = sys.stdin
    if lazy:
        suppl = generate_mols_from_json(StreamJsonListLoader(input))
    else:
        suppl = generate_mols_from_json(json.load(input))
    return input, suppl


def default_open_output(outputDef, defaultOutput, outputFormat):
    if not outputFormat:
        log("No output format specified - using sdf")
        outputFormat = 'sdf'
    if not outputDef:
        outputBase = defaultOutput
    else:
        outputBase = outputDef

    if outputFormat == 'sdf':
        output,writer = default_open_output_sdf(outputDef, outputBase)
    elif outputFormat == 'json':
        output,writer = default_open_output_json(outputDef, outputBase)
    else:
        raise ValueError('Unsupported output format')
    return output,writer,outputBase



def write_basic_squonk_datasetmetadata_hack(outputBase):
    """This is a temp hack to write the minimal metadata that Squonk needs.
    Will needs to be replaced with something that allows something more complete to be written.

    :param outputBase: Base name for the file to write to
    """
    d = {}
    d['type'] = 'org.squonk.types.MoleculeObject'
    s = json.dumps(d)
    meta = open(outputBase + '.metadata', 'w')
    meta.write(s)
    meta.close()


def default_open_output_sdf(outputDef, outputBase):
    if outputDef:
        output = gzip.open(outputDef + '.sdf.gz','w+')
    else:
        output = sys.stdout
    writer = Chem.SDWriter(output)
    return output,writer


def default_open_output_json(outputDef, outputBase):

    # this is a temp hack write some basic metadata that Squonk needs
    write_basic_squonk_datasetmetadata_hack(outputBase)

    if outputDef:
        output = gzip.open(outputDef + '.data.gz','w+')
    else:
        output = sys.stdout
    writer = JsonWriter(output)
    return output,writer



def write_metrics(baseName, values):
    """Write the metrics data

    :param baseName: The base name of the output files. e.g. extensions will be appended to this base name
    :param values dictionary of values to write
    """
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



def create_mol_from_props(molobj):
    """Function to get the RDKit mol from MoleculeObject JSON

    :param molobj: Python dictionary containing the molecule's properties
    """
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


def generate_mols_from_json(input):
    """Create a supplier of RDKit Mol objects from the json

    :param input: file like object containing the json representation of the molecules
    """
    j=0
    #for item in items(input, "item"):
    for item in input:
        j+=1
        mol = create_mol_from_props(item)
        yield mol



class JsonWriter:

    def __init__(self, file):
        self.file = file
        self.file.write('[')
        self.count = 0

    def write(self, mol):
        d = {}
        d['source'] = Chem.MolToMolBlock(mol)
        d['format'] = 'mol'
        props = mol.GetPropsAsDict()
        if 'uuid' in props:
            d['uuid'] = props['uuid']
            del props['uuid']
        else:
            d['uuid'] = str(uuid.uuid4())
        d['values'] = props
        json_str = json.dumps(d)
        if self.count > 0:
            self.file.write(',')
        self.file.write(json_str)
        self.count += 1

    def close(self):
        self.file.write(']')
        self.file.close()

    def flush(self):
        self.file.flush()
