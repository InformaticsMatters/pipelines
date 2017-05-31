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
import subprocess


from pipelines.utils import utils


def main():
    global PDB_PATH,WRITER,THRESHOLD
    parser = argparse.ArgumentParser(description='Open babel converter')
    utils.add_default_io_args(parser)
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('-pdb', '--pdb_file', help="PDB file for converting")
    parser.add_argument('-prot', '--protonate', help="Set protonation at pH 7.4", default=False)
    parser.add_argument('-off', '--outputfomat', choices=["pdbqt","mol2"])

    args = parser.parse_args()

    if args.outputfomat == "pdbqt":
        out_format = "-opdbqt"
    elif args.outputfomat == "mol2":
        out_format = "-omol2"
    else:
        raise ValueError("Not accepted output format")

    base_args = ["obabel", "-ipdb", args.pdb_file, out_format, "-O", args.output]
    if args.protonate:
        base_args.append("-p")

    # Now run the process
    p = subprocess.Popen(base_args)
    p.communicate()

if __name__ == "__main__":
    main()
