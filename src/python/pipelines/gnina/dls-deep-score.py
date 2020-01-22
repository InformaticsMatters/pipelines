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
#

# Create dir containing ligands.sdf and protein.pdb
# Enter docker container like this:
#   docker run -it --rm --gpus all -v $PWD:/root/train/fragalysis_test_files/work:Z informaticsmatters/deep-app-ubuntu-1604:latest bash
#
# Now inside the container run like this:
#   rm -rf scored_ligands.* work/proteins work/ligands work/outputs && python3 work/dls-deep-score.py -l ligands.sdf -p protein.pdb -w work

import argparse, os, re
from openbabel import pybel

types_file_name = 'inputs.types'
predict_file_name = 'predictions.txt'
work_dir = '.'

def write_inputs(protein_file, ligands_file):
    global types_file_name
    global work_dir
    ligands_path = "{0}/ligands".format(work_dir)
    print("Writing ligands to", ligands_path)
    os.mkdir(ligands_path)
    cmd1 = "gninatyper {0} {1}/ligands/ligand".format(ligands_file, work_dir)
    print('CMD:', cmd1)
    os.system(cmd1)
    ligand_gninatypes = os.listdir("{0}/ligands".format(work_dir))

    proteins_path = "{0}/proteins".format(work_dir)
    print("Writing proteins to", proteins_path)
    os.mkdir(proteins_path)
    cmd2 = "gninatyper {0} {1}/proteins/protein".format(protein_file, work_dir)
    print('CMD:', cmd2)
    os.system(cmd2)
    protein_gninatypes = os.listdir("{0}/proteins".format(work_dir))

    types_path = "{0}/{1}".format(work_dir, types_file_name)
    print("Writing types to", types_path)
    with open(types_path, 'w') as types_file:
        for protein in protein_gninatypes:
            for ligand in ligand_gninatypes:
                line = "0 {0}/proteins/{1} {0}/ligands/{2}\n".format(work_dir, protein, ligand)
                types_file.write(line)

def run_predictions():
    global types_file_name
    global predict_file_name
    global work_dir
    # python3 scripts/predict.py -m resources/dense.prototxt -w resources/weights.caffemodel -i work_0/test_set.types >> work_0/caffe_output/predictions.txt
    cmd1 = "python3 /root/train/fragalysis_test_files/scripts/predict.py -m /root/train/fragalysis_test_files/resources/dense.prototxt" +\
           " -w /root/train/fragalysis_test_files/resources/weights.caffemodel" +\
            " -i {0}/{1} -o {0}/{2}".format(work_dir, types_file_name, predict_file_name)
    print("CMD:", cmd1)
    os.system(cmd1)

def read_predictions(ligands_file):
    global predict_file_name
    global work_dir
    scores = {}
    with open("{0}/{1}".format(work_dir, predict_file_name), 'r') as input:
        for line in input:
            #print(line)
            tokens = line.split()
            if len(tokens) == 5 and tokens[1] == '|':
                # print(len(tokens), tokens[0], tokens[3], tokens[4])
                record_no = match_ligand(tokens[4])
                if record_no is not None:
                    # print(record_no, tokens[0])
                    scores[record_no] = tokens[0]
    print("Found", len(scores), "scores")
    return scores

def patch_scores(sdf_in, outfile, scores):

    global work_dir

    counter = 0
    sdf_path = "{0}/{1}.sdf".format(work_dir, outfile)
    tsv_path = "{0}/{1}.tsv".format(work_dir, outfile)
    print("Writing results to {0} and {1}".format(sdf_path, tsv_path))
    with open(tsv_path, 'w') as tsv_file:
        sdf_file = pybel.Outputfile("sdf", sdf_path)
        for mol in pybel.readfile("sdf", sdf_in):
            if counter in scores:
                score = scores[counter]
                # print("Score for record {0} is {1}".format(counter, score))

                mol.data['dls_deep_score'] = score
                if 'SCORE' in mol.data:
                    rdock_score = mol.data['SCORE']
                else:
                    rdock_score = ''

                if 'SCORE.norm' in mol.data:
                    rdock_nscore = mol.data['SCORE.norm']
                else:
                    rdock_nscore = ''

                sdf_file.write(mol)
                tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(counter, rdock_score, rdock_nscore, score))

            else:
                print("No score found for record", counter)
            counter += 1
        sdf_file.close()

# work/ligands/ligand_9.gninatypes
ligand_patt = re.compile(r'.*ligands/ligand_(\d+)\.gninatypes$')

def match_ligand(s):
    global ligand_patt
    m = ligand_patt.match(s)
    if m:
        i = m.group(1)
        return int(i)
    else:
        return None

def main():

    global work_dir


    parser = argparse.ArgumentParser(description='DLS Deep - pose scoring')
    # parameter_utils.add_default_io_args(parser)
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('-l', '--ligands', help="Ligands to be scored")
    parser.add_argument('-p', '--pdb-file', help="PDB file for scoring")
    parser.add_argument('-o', '--outfile', default='scored_ligands', help="Base file name for results")
    parser.add_argument('-w', '--work-dir', default=".", help="Working directory")
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()
    print("DLS deep args: ", args)

    work_dir = args.work_dir

    # if args.informat != 'sdf':
        # TODO convert to SDF
    # # Open up the input file
    # input, suppl = rdkit_utils.default_open_input(args.input, args.informat)
    # # Open the output file
    # s_now = datetime.datetime.utcnow().strftime("%d-%b-%Y %H:%M:%S UTC")
    # source = 'pipelines/gnina/dls-deep-score.py'
    # output, WRITER, output_base = \
    #     rdkit_utils.default_open_output(args.output, "dls-deep-score", args.outformat,
    #                                     compress=not args.no_gzip,
    #                                     thinOutput=args.thin,
    #                                     valueClassMappings={'dls-deep-score': 'java.lang.Float'},
    #                                     datasetMetaProps={'created': s_now,
    #                                                       'source': source,
    #                                                       'description': 'DLS Deep - pose scoring'}
    #                                     )
    #
    # PDB_PATH = args.pdb_file
    #    # Close the file
    # WRITER.close()

    protein_pdb = args.pdb_file
    ligands_sdf = args.ligands
    outfile = args.outfile
    write_inputs(protein_pdb, ligands_sdf)
    run_predictions()
    scores = read_predictions(ligands_sdf)
    patch_scores(ligands_sdf, outfile, scores)

    # if args.meta:
    #     utils.write_metrics(output_base, {'__InputCount__': COUNTER, '__OutputCount__': SUCCESS, 'PLI': COUNTER})


if __name__ == "__main__":
    main()
