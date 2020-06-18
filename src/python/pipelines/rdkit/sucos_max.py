#!/usr/bin/env python
"""
Assess ligands against a second set of molecules using SuCOS scores.
This is a quite specialised function that is designed to take a set of potential follow up
compounds and compare them to a set of clustered fragment hits to help identify which follow up
ligands best map to the binding space of the hits.

The clustering of the fragment hits is expected to be performed with the sucos_cluster.py module
and will generate a set of SD files, one for each cluster of hits (presumably corresponding to a
binding pocket in the protein target).

Each molecule in the input ligands is then compared (using SuCOS) to each hit in the clusters. There
 are different modes which determine how the ligand is assessed.

In mode 'max' the hit with the best SuCOS score is identified. The output is a SD file with each of the ligands,
with these additional fields for each molecule:
Max_SuCOS_Score - the best score
Max_SuCOS_FeatureMap_Score - the feature map score for the hit that has the best SuCOS score
Max_SuCOS_Protrude_Score - the protrude volume for the hit that has the best SuCOS score
Max_SuCOS_Cluster - the name of the cluster SD file that contains the best hit
Max_SuCOS_Index - the index of the best hit in the SD file

In mode 'cum' the sum of all the scores is calculated and reported as the following properties for each molecule:
Cum_SuCOS_Score property: the sum of the SuCOS scores
Cum_SuCOS_FeatureMap_Score: the sum of the feature map scores
Cum_SuCOS_Protrude_Score: the sum of the protrude volume scores

If a molecule has no alignment to any of the clustered hits (all alignment scores of zero) then it is not
included in the results.


SuCOS is the work of Susan Leung.
GitHub: https://github.com/susanhleung/SuCOS
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""

import sys, argparse, gzip, os
from rdkit import Chem, rdBase
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils
from . import sucos


field_SuCOSMax_Score = "SuCOS_Max_Score"
field_SuCOSMax_FMScore = "SuCOS_Max_FeatureMap_Score"
field_SuCOSMax_ProtrudeScore = "SuCOS_Max_Protrude_Score"
field_SuCOSMax_Index = "SuCOS_Max_Index"
field_SuCOSMax_Target = "SuCOS_Max_Target"

field_SuCOSCum_Score = "SuCOS_Cum_Score"
field_SuCOSCum_FMScore = "SuCOS_Cum_FeatureMap_Score"
field_SuCOSCum_ProtrudeScore = "SuCOS_Cum_Protrude_Score"


def process(inputs_supplr, targets_supplr, writer, field_name, filter_value, filter_field):
    cluster = []
    mol_ids = []
    i = 0
    for mol in targets_supplr:
        i += 1
        if not mol:
            utils.log("WARNING: failed to generate target molecule", i)
            continue
        try:
            if field_name:
                id = mol.GetProp(field_name)
                mol_ids.append(id)
            else:
                mol_ids.append(None)
            features = sucos.getRawFeatures(mol)
            cluster.append((mol, features))
        except:
            utils.log("WARNING: failed to generate features for molecule", i, sys.exc_info())
    utils.log("Generated features for", len(cluster), "molecules")

    comparisons = 0
    mol_num = 0
    errors = 0

    for mol in inputs_supplr:
        mol_num += 1
        if not mol:
            utils.log("WARNING: failed to generate molecule", mol_num, "in input")
            errors += 1
            continue
        try:
            query_features = sucos.getRawFeatures(mol)
        except:
            utils.log("WARNING: failed to generate features for molecule", mol_num, "in input")
            errors += 1
            continue

        scores_max = [0, 0, 0]
        scores_cum = [0, 0, 0]

        index = 0
        for entry in cluster:
            hit = entry[0]
            ref_features = entry[1]
            comparisons += 1
            sucos_score, fm_score, vol_score = sucos.get_SucosScore(hit, mol,
                                                                        tani=False, ref_features=ref_features,
                                                                        query_features=query_features)

            if sucos_score > scores_max[0]:
                scores_max[0] = sucos_score
                scores_max[1] = fm_score
                scores_max[2] = vol_score
                cluster_index = index
                best_id = mol_ids[index]

            scores_cum[0] += sucos_score
            scores_cum[1] += fm_score
            scores_cum[2] += vol_score

            index += 1

        # utils.log("Max SuCOS:", scores[0], "FM:", scores[1], "P:", scores[2],"File:", cluster_file_name_only, "Index:", cluster_index)
        mol.SetDoubleProp(field_SuCOSMax_Score, scores_max[0] if scores_max[0] > 0 else 0)
        mol.SetDoubleProp(field_SuCOSMax_FMScore, scores_max[1] if scores_max[1] > 0 else 0)
        mol.SetDoubleProp(field_SuCOSMax_ProtrudeScore, scores_max[2] if scores_max[2] > 0 else 0)

        if best_id:
            mol.SetProp(field_SuCOSMax_Target, best_id)
            mol.SetIntProp(field_SuCOSMax_Index, cluster_index)

        # utils.log("Cum SuCOS:", scores[0], "FM:", scores[1], "P:", scores[2])
        mol.SetDoubleProp(field_SuCOSCum_Score, scores_cum[0] if scores_cum[0] > 0 else 0)
        mol.SetDoubleProp(field_SuCOSCum_FMScore, scores_cum[1] if scores_cum[1] > 0 else 0)
        mol.SetDoubleProp(field_SuCOSCum_ProtrudeScore, scores_cum[2] if scores_cum[2] > 0 else 0)


        if filter_value and filter_field:
            if mol.HasProp(filter_field):
                val = mol.GetDoubleProp(filter_field)
                if val > filter_value:
                    writer.write(mol)
        else:
            writer.write(mol)

    utils.log("Completed", comparisons, "comparisons")
    return mol_num, comparisons, errors


### start main execution #########################################

def main():
    parser = argparse.ArgumentParser(description='Max SuCOS scores with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-tm', '--target-molecules', help='Target molecules to compare against')
    parser.add_argument('-tf', '--targets-format', help='Target molecules format')
    parser.add_argument('-n', '--name-field', help='Name of field with molecule name')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('--filter-value', type=float, help='Filter out values with scores less than this.')
    parser.add_argument('--filter-field', help='Field to use to filter values.')

    args = parser.parse_args()
    utils.log("Max SuCOSMax Args: ", args)

    source = "sucos_max.py"
    datasetMetaProps = {"source": source, "description": "SuCOSMax using RDKit " + rdBase.rdkitVersion}
    clsMappings = {}
    fieldMetaProps = []

    clsMappings[field_SuCOSMax_Score] = "java.lang.Float"
    clsMappings[field_SuCOSMax_FMScore] = "java.lang.Float"
    clsMappings[field_SuCOSMax_ProtrudeScore] = "java.lang.Float"
    clsMappings[field_SuCOSMax_Index] = "java.lang.Integer"
    clsMappings[field_SuCOSCum_Score] = "java.lang.Float"
    clsMappings[field_SuCOSCum_FMScore] = "java.lang.Float"
    clsMappings[field_SuCOSCum_ProtrudeScore] = "java.lang.Float"

    fieldMetaProps.append(
        {"fieldName": field_SuCOSMax_Score, "values": {"source": source, "description": "SuCOS Max score"}})
    fieldMetaProps.append({"fieldName": field_SuCOSMax_FMScore,
                           "values": {"source": source, "description": "SuCOS Max Feature Map score"}})
    fieldMetaProps.append({"fieldName": field_SuCOSMax_ProtrudeScore,
                           "values": {"source": source, "description": "SuCOS Max Protrude score"}})
    fieldMetaProps.append(
        {"fieldName": field_SuCOSMax_Index, "values": {"source": source, "description": "SuCOS Max target index"}})
    fieldMetaProps.append(
        {"fieldName": field_SuCOSCum_Score, "values": {"source": source, "description": "SuCOS Cumulative score"}})
    fieldMetaProps.append({"fieldName": field_SuCOSCum_FMScore,
                           "values": {"source": source, "description": "SuCOS Cumulative Feature Map score"}})
    fieldMetaProps.append({"fieldName": field_SuCOSCum_ProtrudeScore,
                           "values": {"source": source, "description": "SuCOS Cumulative Protrude score"}})

    if args.name_field:
        clsMappings[field_SuCOSMax_Target] = "java.lang.String"
        fieldMetaProps.append(
            {"fieldName": field_SuCOSMax_Target, "values": {"source": source, "description": "SuCOS Max target name"}})

    inputs_file, inputs_supplr = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                                                                  'sucos-max', args.outformat,
                                                                  valueClassMappings=clsMappings,
                                                                  datasetMetaProps=datasetMetaProps,
                                                                  fieldMetaProps=fieldMetaProps,
                                                                  compress=not args.no_gzip)

    targets_file, targets_supplr = rdkit_utils.default_open_input(args.target_molecules, args.targets_format)

    count, total, errors = process(inputs_supplr, targets_supplr, writer, args.name_field, args.filter_value, args.filter_field)

    inputs_file.close()
    targets_file.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': count, '__OutputCount__': total, '__ErrorCount__': errors,
                                          'RDKitSuCOS': total})


if __name__ == "__main__":
    main()
