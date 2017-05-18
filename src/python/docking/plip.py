#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
import subprocess

from rdkit import Chem

from src.python import utils


def run_and_get_ans(mol, pdb_path):
    smogmol = '/tmp/smogmol.sdf'
    out_f = open(smogmol, "w")
    out_f.write(Chem.MolToMolBlock(mol))
    out_f.close()
    # Run command
    pli_path = "/usr/local/pli/bin/pli"
    cmd= [pli_path,"-protein",pdb_path,"-ligand",smogmol,"-mode","score","-output","system,scores","-exact_voronoi_areas","0","-selection","ligand","-oformat","json","-minimise","1","-warnings","0","-min_max_iter", "10"]
    proc = subprocess.Popen(cmd,stdout=subprocess.PIPE)
    # Parse the output
    me = proc.stdout.read()
    if not me:
        return None
    return json.loads(me)

def main():

    parser = argparse.ArgumentParser(description='SMoG2016 - Docking calculation.')
    utils.add_default_io_args(parser)
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('-pdb', '--pdb_file', help="PDB file for scoring")
    args = parser.parse_args()

    # Open up the input file
    input, suppl = utils.default_open_input(args.input, args.informat)
    # Open the ouput file

    output, writer, output_base = utils.default_open_output(args.output, "plip", args.outformat, compress=not args.no_gzip)

    # Iterate over the molecules
    for mol in suppl:
        answer_dict = run_and_get_ans(mol, args.pdb_file)
        if not answer_dict:
            utils.log("FAILED MOL", Chem.MolToSmiles(mol))
            continue
        for ans in answer_dict["system"]:
            if ans.startswith(u"pliff"):
                mol.SetDoubleProp(str(ans),answer_dict["system"][ans])
        utils.log("SCORED MOL:", Chem.MolToSmiles(mol), answer_dict)
        # Write ligand
        writer.write(mol)
        writer.flush()
    # Close the file
    writer.close()

if __name__ == "__main__":
    main()
