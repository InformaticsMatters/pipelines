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
import os,shutil
import subprocess
import tempfile
import threading
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool

from rdkit import Chem

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

lock = threading.Lock()
PDB_PATH = ""
WRITER = ""
COUNTER = 0
SUCCESS = 0

def run_and_get_ans(mol):
    global PDB_PATH
    smogmol = tempfile.NamedTemporaryFile("w",suffix=".sdf",delete=True).name
    out_f = open(smogmol, "w")
    out_f.write(Chem.MolToMolBlock(mol))
    out_f.close()
    # Run command
    proc = subprocess.Popen(["/usr/local/SMoG2016_Rev1/SMoG2016.exe", PDB_PATH, smogmol, "DeltaG"],
                            stdout=subprocess.PIPE)
    # Parse the output
    me = proc.stdout.read()
    if not me:
        # TODO - shouldn't we fail instead?
        return None
    answer = float(me.split("DeltaG")[1].strip())
    return answer

def run_dock(mol):
    global COUNTER
    global SUCCESS
    global THRESHOLD
    answer = run_and_get_ans(mol)
    COUNTER+=1
    if answer is None:
        utils.log("FAILED MOL", Chem.MolToSmiles(mol))
        return
    if THRESHOLD is not None:
        print(answer,THRESHOLD)
        if answer > THRESHOLD:
            utils.log("UNDER THRESHOLD", Chem.MolToSmiles(mol))
            return
    mol.SetDoubleProp("SMoG2016_SCORE", answer)
    utils.log("SCORED MOL:", Chem.MolToSmiles(mol), answer)
    # Write ligand
    lock.acquire()
    SUCCESS+=1
    WRITER.write(mol)
    WRITER.flush()
    lock.release()
    return

def main():
    global WRITER,THRESHOLD
    global PDB_PATH
    parser = argparse.ArgumentParser(description='SMoG2016 - Docking calculation.')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('-pdb', '--pdb_file', help="PDB file for scoring")
    parser.add_argument('-t', '--threshold', help="The maximum score to allow", default=None)
    parser.add_argument('--threads', type=int, help="Number of threads to used. Default is the number of cores", default=None)
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()

    utils.log("SMoG2016 Args: ", args)

    smog_path = "/usr/local/SMoG2016_Rev1/"
    if args.threshold:
        THRESHOLD = float(args.threshold)
    else:
        THRESHOLD = None

    PDB_PATH = "/tmp/pdb_file.pdb"
    # Now copy it to prot_pdb.pdb -> silly SMOG bug requires underscore in the filename!
    shutil.copy(args.pdb_file,PDB_PATH)

    # Open up the input file
    input, suppl = rdkit_utils.default_open_input(args.input, args.informat)
    # Open the output file
    output, WRITER, output_base = rdkit_utils.\
        default_open_output(args.output, "SMoG2016",
                            args.outformat, compress=not args.no_gzip)

    # Cd to the route of the action
    # TODO - can this be done without changing dir? It gives problems in finding the input files and in writing the metrics
    cwd = os.getcwd()
    os.chdir(smog_path)

    # Iterate over the molecules
    # WARNING - if using parallel processing the order of molecules is not preserved. Set args.threads to 1 to ensure this.
    if args.threads is None:
        threads = multiprocessing.cpu_count()
    else:
        threads = args.threads
    pool = ThreadPool(threads)
    pool.map(run_dock, suppl)
    # Close the file
    WRITER.close()

    os.chdir(cwd)
    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': COUNTER, '__OutputCount__': SUCCESS, 'SMoG2016': COUNTER})

    utils.log("SMoG2016 complete")

if __name__ == "__main__":
    main()
