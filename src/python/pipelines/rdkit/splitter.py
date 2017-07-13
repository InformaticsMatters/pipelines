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

from pipelines.utils import utils


def main():

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='RDKit Input Splitter')
    utils.add_default_input_args(parser)
    parser.add_argument('-o', '--output', required=True, help="Directory name for output files (no extension).")
    parser.add_argument('-f', '--field', required=True, help="field to use to split input. Output files will have the name of this field's value")
    parser.add_argument('--meta', action='store_true', help='Write metadata and metrics files')

    args = parser.parse_args()
    utils.log("Splitter Args: ", args)

    filenames = split(args.input, args.informat, args.field, args.output, args.meta)
    utils.log("Files generated:", " ".join(filenames))


def split(input, informat, fieldName, outputBase, writeMetrics):
    """Splits the input into separate files. The name of each file and the file the each record is written to
    is determined by the fieldName parameter
    """

    input,suppl = utils.default_open_input(input, informat)

    i=0
    written=0
    writers = {}
    outputs = []
    filenames = []
    for mol in suppl:
        i +=1
        if mol is None: continue
        if not mol.HasProp(fieldName):
            utils.log("Skipping molecule", i, "- did not contain field",fieldName)
            continue
        value = mol.GetProp(fieldName)
        if value:
            s = str(value)
            if writers.has_key(s):
                writer = writers[s]
            else:
                name = outputBase + s
                output, writer = utils.default_open_output_sdf(name, outputBase, False, False)
                filenames.append(name + '.sdf')
                outputs.append(output)
                writers[s] = writer
            writer.write(mol)
            written += 1


    utils.log("Generated", len(writers), "outputs from", i, "records")

    input.close()
    for k in writers: writers[k].close()
    for o in outputs: o.close()

    if writeMetrics:
        utils.write_metrics(outputBase, {'__InputCount__':i, '__OutputCount__':written, 'Splitter':i})

    return filenames

if __name__ == "__main__":
    main()