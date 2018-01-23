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
import sys, subprocess

from pipelines_utils import utils


def execute(input, output, extension, format, ph, noGzip):
    filename = output + "." + extension
    base_args = ["obabel", "-ipdb", input, format, "-O", filename]
    if ph:
        base_args.append("-p")
        base_args.append("ph")

    subprocess.check_call(base_args, stdout=sys.stderr, stderr=sys.stderr)

    # NOTE the -z argument does not seem to work correctly with obabel (truncted files generated) so we
    # fall back to good old gzip to handle the compression once the uncompressed file is created
    if not noGzip:
        subprocess.check_call(['gzip', filename], stdout=sys.stderr, stderr=sys.stderr)

def main():
    global PDB_PATH,WRITER,THRESHOLD
    parser = argparse.ArgumentParser(description='Open babel PDB prepare')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output')
    parser.add_argument('-i', '--input', help="PDB file for converting")
    parser.add_argument('-o', '--output', help="Base name for output files (no extension).")
    parser.add_argument('-mol2', '--mol2', action='store_true', help='Output as Mol2 format.')
    parser.add_argument('-pdbqt', '--pdbqt', action='store_true', help='Output as pdbqt format.')
    parser.add_argument('--meta', action='store_true', help='Write metrics files')
    parser.add_argument('-prot', '--protonate', type=float, help="protonate at this pH (optional)")

    args = parser.parse_args()

    utils.log("Prepare Args: ", args)

    if not (args.mol2 or args.pdbqt):
        raise ValueError("Must specify at least one output fromat: mol2 and/or pdbqt")


    if args.pdbqt:
        utils.log("Preparing as pdbqt")
        execute(args.input, args.output, "pdbqt", "-opdbqt", args.protonate, args.no_gzip)

    if args.mol2:
        utils.log("Preparing as mol2")
        execute(args.input, args.output, "mol2", "-omol2", args.protonate, args.no_gzip)

    utils.log("Preparation complete")


if __name__ == "__main__":
    main()
