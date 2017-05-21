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

    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit Sdf2Json')
    parser.add_argument('-i', '--input', help="Input SD file, if not defined the STDIN is used")
    parser.add_argument('-o', '--output', help="Base name for output json file (no extension). If not defined then SDTOUT is used for the structures and output is used as base name of the other files.")
    parser.add_argument('--exclude', help="Optional list of fields (comma separated) to exclude from the output.")


    args = parser.parse_args()
    utils.log("Screen Args: ", args)

    if args.input:
        if args.input.lower().endswith(".sdf"):
            base = args.input[:-4]
        elif args.input.lower().endswith(".sdf.gz"):
            base = args.input[:-7]
        else:
            base = "json"
    utils.log("Base:", base)


    input,output,suppl,writer,output_base = utils.default_open_input_output(args.input, "sdf", args.output, base, "json")
    if args.exclude:
        excludes = args.exclude.split(",")
        utils.log("Excluding", excludes)
    else:
        excludes = None

    i=0
    count = 0
    for mol in suppl:
        i +=1
        if mol is None: continue
        if excludes:
            for exclude in excludes:
                if mol.HasProp(exclude): mol.ClearProp(exclude)
        writer.write(mol)
        count += 1

    utils.log("Converted", count, " molecules")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    utils.write_metrics(output_base, {'__InputCount__':i, '__OutputCount__':count, 'RDKitSdf2Json':count})

    return count
    
if __name__ == "__main__":
    main()

