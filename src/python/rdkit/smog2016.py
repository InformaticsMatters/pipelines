import os
from rdkit import Chem
import utils,argparse,subprocess


def run_and_get_ans(mol, mol_path, pdb_path):
    out_f = open(mol_path, "w")
    out_f.write(Chem.MolToMolBlock(mol))
    out_f.close()
    # Run command
    proc = subprocess.Popen(["/usr/local/SMoG2016_Rev1/SMoG2016.exe", pdb_path, mol_path, "DeltaG"],
                            stdout=subprocess.PIPE)
    # Parse the output
    me = proc.stdout.read()
    if not me:
        return 100.0
    answer = float(me.split("DeltaG")[1].strip())
    return answer

def main():


    parser = argparse.ArgumentParser(description='SMoG2016 - Docking calculation.')
    utils.add_default_io_args(parser)
    parser.add_argument('-pdb', '--pdb_file', help="PDB file for scoring")
    args = parser.parse_args()

    smog_path = "/usr/local/SMoG2016_Rev1/"
    # Cd to the route of the action
    os.chdir(smog_path)

    # Open up the input file
    input, suppl = utils.default_open_input(args.input, args.informat)
    # Open the ouput file
    output, writer, output_base = utils.default_open_output(args.output, "SMoG2016", args.outformat)
    out_sd = Chem.SDWriter(args.output)

    input_path = os.path.join(smog_path,"mol.sdf")
    # Iterate over the molecules
    for mol in suppl:
        print("SCORING MOL")
        answer = run_and_get_ans(mol,input_path, args.pdb_file)
        mol.SetDoubleProp("SMoG2016_SCORE",answer)
        print("SCORED MOL:"+str(answer))
        print("SCORED MOL:"+Chem.MolToSmiles(mol))
        # Write ligand
        writer.write(mol)
        out_sd.write(mol)
        writer.flush()
    # Close the file
    writer.close()

if __name__ == "__main__":
    main()
