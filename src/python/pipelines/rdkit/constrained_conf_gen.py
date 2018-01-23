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

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.MCS import FindMCS

from pipelines_utils import utils
import sys, logging


def guess_substruct(mol_one, mol_two):
    """Code to find the substructure between two molecules."""
    return Chem.MolToSmiles(Chem.MolFromSmarts(FindMCS([mol_one,mol_two],completeRingsOnly=True,matchValences=True).smarts))

def generate_conformers(molIdx, my_mol, NumOfConf, ref_mol, outputfile, coreSubstruct):
    # Find the MCS if not given
    if not coreSubstruct:
        coreSubstruct = guess_substruct(my_mol,ref_mol)

    # Creating core of reference ligand #
    core1 = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(ref_mol, Chem.MolFromSmiles(coreSubstruct)),
                                     Chem.MolFromSmiles('*'))
    core1.UpdatePropertyCache()

    # Generate conformers with constrained embed
    conf_lst = []
    count = 0
    errors = 0
    for i in (xrange(NumOfConf)):
        try:
            conf_lst.append(Chem.AddHs(my_mol))
            AllChem.ConstrainedEmbed(conf_lst[i], core1, randomseed=i)
            cleaned = Chem.RemoveHs(conf_lst[i])
            cleaned.ClearProp("uuid")
            if my_mol.HasProp("uuid"):
                cleaned.SetProp("SourceMolUUID", my_mol.GetProp("uuid"))
            cleaned.SetIntProp("SourceMolNum", molIdx)
            cleaned.SetIntProp("ConformerNum", count + 1)
            outputfile.write(cleaned)
            count+=1
        except ValueError:
            errors+=1
            logging.exception('')
    return count, errors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RDKit constrained conformer generator')
    utils.add_default_io_args(parser)
    parser.add_argument('-n', '--num', type=int, default=10, help='number of conformers to generate')
    parser.add_argument('-r', '--refmol', help="Reference molecule file")
    parser.add_argument('--refmolidx', help="Reference molecule index in file", type=int, default=1)
    parser.add_argument('-c', '--core_smi', help='Core substructure. If not specified - guessed using MCS', default='')


    args = parser.parse_args()
    # Get the reference molecule
    ref_mol_input, ref_mol_suppl = utils.default_open_input(args.refmol, args.refmol)
    counter = 0
    # Get the specified reference molecule. Default is the first
    for mol in ref_mol_suppl:
        counter+=1
        if counter == args.refmolidx:
            ref_mol = mol
            break
    ref_mol_input.close()

    if counter < args.refmolidx:
        raise ValueError("Invalid refmolidx. " + str(args.refmolidx) + " was specified but only " + str(counter) + " molecules were present in refmol.")


    # handle metadata
    source = "constrained_conf_gen.py"
    datasetMetaProps = {"source":source, "description": "Constrained conformer generation using RDKit " + rdBase.rdkitVersion}
    clsMappings = {"EmbedRMS": "java.lang.Float"}
    fieldMetaProps = [{"fieldName":"EmbedRMS", "values": {"source":source, "description":"Embedding RMS value"}}]

    # Get the molecules
    input, suppl = utils.default_open_input(args.input, args.informat)
    output, WRITER, output_base = utils.default_open_output(args.output, "const_conf_gen", args.outformat,
                                                            valueClassMappings=clsMappings, datasetMetaProps=datasetMetaProps, fieldMetaProps=fieldMetaProps)

    inputs = 0
    totalCount = 0
    totalErrors = 0
    for mol in suppl:
        inputs += 1
        if mol:
            count, errors = generate_conformers(inputs, mol, args.num, ref_mol, WRITER, args.core_smi)
            totalCount += count
            totalErrors += errors

    input.close()
    WRITER.close()

    if totalErrors > 0:
        utils.log("WARNING:", totalErrors, "conformers failed to generate")

    # write metrics
    if args.meta:
        metrics = {'__InputCount__':inputs, '__OutputCount__':totalCount, 'RDKitConstrainedConformer':totalCount}
        if totalErrors > 0:
            metrics['__ErrorCount__'] = totalErrors
        utils.write_metrics(output_base, metrics)
