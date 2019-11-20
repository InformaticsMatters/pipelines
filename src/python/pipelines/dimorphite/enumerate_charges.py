import argparse, sys

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils

from .dimorphite_dl import run_with_mol_list

def enumerateMol(mol, fragment):
    """
    Enumerate a single molecule
    :param mol:
    :param fragment The fragmentation method, 'hac' or 'mw'. If not specified the whole molecules is passed to Dimorphite
    :return:
    """

    if fragment:
        mol = mol_utils.fragment(mol, fragment)

    inputmol = []
    inputmol.append(mol)

    protonated_mols = run_with_mol_list(inputmol)
    return protonated_mols


def writeEnumeratedMols(src_mol, enum_mols, writer, index):
    """
    Write the enumerated molecule to the writer, incorporating the ID of the source molecule.
    :param src_mol:
    :param enum_mols:
    :param writer:
    :param index:
    :return:
    """
    errors = 0
    total = 0
    for mol in enum_mols:
        if mol is None:
            errors += 1
            continue
        add_src_mol_ref(src_mol, mol, index)
        writer.write(mol)
        total += 1

    return total, errors


def add_src_mol_ref(src_mol, target_mol, index):
    """
    Add the ID of the source molecule to the enumerated molecule as the field named EnumChargeSrcMol.
    The ID is taken form the uuid field if it exists, if not form the _Name field if it exists and finally
    from the index parameter (the index of the source molecule in the input) if neither of those fields are found.
    :param src_mol:
    :param target_mol:
    :param index:
    :return:
    """
    if src_mol.HasProp('uuid'):
        parent = src_mol.GetProp('uuid')
    elif src_mol.HasProp('_name_'):
        parent = src_mol.GetProp('_Name')
    else:
        parent = str(index)

    if parent:
        target_mol.SetProp('EnumChargeSrcMol', parent)

### start main execution #########################################

def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Enumerate charges')
    parser.add_argument('--fragment-method', choices=['hac', 'mw'],
                        help='Approach to find biggest fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight)')
    parser.add_argument('--min-ph', help='The min pH to consider', type=float, default=5.0)
    parser.add_argument('--max-ph', help='The max pH to consider', type=float, default=9.0)

    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()
    utils.log("Enumerate charges: ", args)

    # handle metadata
    source = "enumerate_charges.py"
    datasetMetaProps = {"source":source, "description": "Enumerate charges using Dimorphite-dl"}
    clsMappings = {
        "EnumChargeSrcMol": "java.lang.String"}
    fieldMetaProps = [
        {"fieldName":"EnumChargeSrcMol",   "values": {"source":source, "description":"ID of source molecule"}}
    ]

    oformat = utils.determine_output_format(args.outformat)

    input,output,suppl,writer,output_base = rdkit_utils. \
        default_open_input_output(args.input, args.informat, args.output,
                                  'enumerateCharges', args.outformat,
                                  thinOutput=False, valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)

    count = 0
    total = 0
    errors = 0
    min_ph = args.min_ph
    max_ph = args.max_ph

    # this hacky bit is needed because the dimporphite entrypoint assumes it's args are passed using argparse
    # but it doesn't understand our args, so we need to switch between the two sets of args.
    dimorphite_sys_argv = sys.argv[:1]
    dimorphite_sys_argv.append('--min_ph')
    dimorphite_sys_argv.append(str(min_ph))
    dimorphite_sys_argv.append('--max_ph')
    dimorphite_sys_argv.append(str(max_ph))
    fragment = args.fragment_method
    for mol in suppl:
        count += 1
        orig_sys_argv = sys.argv[:]
        sys.argv = dimorphite_sys_argv
        enum_mols = enumerateMol(mol, fragment)
        sys.argv = orig_sys_argv
        t, e = writeEnumeratedMols(mol, enum_mols, writer, count)
        total += t
        errors += e

    utils.log(count, total, errors)

    input.close()
    writer.flush()
    writer.close()
    output.close()

    # re-write the metadata as we now know the size
    if oformat == 'json':
        utils.write_squonk_datasetmetadata(output_base, False, clsMappings, datasetMetaProps, fieldMetaProps, size=total)

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':count, '__OutputCount__':total, '__ErrorCount__':errors, 'EnumerateChargesDimporphite':total})

if __name__ == "__main__":
    main()
