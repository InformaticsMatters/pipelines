import argparse, sys

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils, mol_utils
from rdkit import Chem

from . import calc_interactions


def execute(suppl, writer, report, group_by_field, score_field, score_descending, stats_fields):
    count = 0
    total = 0
    errors = 0
    group_count = 0
    curr_gbv = None
    interactions = {}
    best_scores = {}
    best_mols =  {}
    stats_data = {}
    for mol in suppl:
        total +=1
        if not mol:
            errors += 1
            continue
        if not mol.HasProp(group_by_field):
            report.write("WARNING: molecule %s does not contain field %s\n" % (total, group_by_field))
            errors += 1
            continue
        if not mol.HasProp(score_field):
            report.write("WARNING: molecule %s does not contain field %s\n" % (total, score_field))
            errors += 1
            continue
        gbv = mol.GetProp(group_by_field)
        sco = mol.GetDoubleProp(score_field)
        inters = gen_interactions(mol)
        utils.log('processing', gbv, inters)
        if gbv != curr_gbv:

            # write summaries
            if curr_gbv:
                write_summary(writer, report, gbv, group_count, best_mols, stats_data)
            curr_gbv = gbv
            group_count = 1
            best_scores = {inters: sco}
            best_mols = {inters: mol}
            stats_data = {}
            add_stats(mol, inters, stats_fields, stats_data)
        else:
            # add to summary
            group_count += 1
            curr_best_sco = best_scores.get(inters, None)
            add_stats(mol, inters, stats_fields, stats_data)
            if None == curr_best_sco:
                best_scores[inters] = sco
                best_mols[inters] = mol
            else:
                if score_descending:
                    if sco > curr_best_sco:
                        best_scores[inters] = sco
                        best_mols[inters] = mol
                else:
                    if sco < curr_best_sco:
                        best_scores[inters] = sco
                        best_mols[inters] = mol
        count += 1
    write_summary(writer, report, gbv, group_count, best_mols, stats_data)

    return count, total, errors

def write_summary(writer, report, gbv, count, best_mols, stats_data):
    report.write("Summary for %s molecules for %s\n" % (count, gbv))
    for inter, mol in best_mols.items():
        report.write("  %s\n" % (str(inter)))
        if inter in stats_data:
            for field, values in stats_data[inter].items():
                report.write("    %s = [%s, %s, %s, %s]\n" % (field, len(values), min(values), max(values), sum(values) / len(values)))
            writer.write(mol)

str_interactions = 'Interactions'

def gen_interactions(mol):
    interactions = []
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_hbond + str_interactions))
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_halogen + str_interactions))
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_hydrophobic + str_interactions))
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_salt_bridge + str_interactions))
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_pi_stacking + str_interactions))
    interactions.append(find_canonical_interactions(mol, calc_interactions.inter_type_pi_cation + str_interactions))
    return tuple(interactions)

def find_canonical_interactions(mol, prop):
    if mol.HasProp(prop):
        canons = []
        inters = mol.GetProp(prop)
        lines = inters.split('\n')
        for line in lines:
            tokens = line.split(' ')
            canon = tokens[0]
            canons.append(canon)
        return tuple(sorted(canons))
    else:
        return None

def add_stats(mol, inters, stats_fields, stats_data):
    if inters in stats_data:
        d = stats_data[inters]
    else:
        d = {}
        stats_data[inters] = d
    if stats_fields:
        for field in stats_fields:
            if mol.HasProp(field):
                v = mol.GetDoubleProp(field)
                if field in d:
                    d[field].append(v)
                else:
                    d[field] = [v]

### start main execution #########################################


def main():
    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Filter interactions')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-f', '--group-by-field', required=True, help='Field to group records by (must be sequential)')
    parser.add_argument('-s', '--score-field', required=True, help='Field to use to rank records within a group')
    parser.add_argument('-d', '--score-descending', action='store_true', help='Sort records in descending order')
    parser.add_argument('-x', '--stats-fields', nargs='*', help='Field to use to for summary statistics')

    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')

    args = parser.parse_args()
    utils.log("filter_interactions: ", args)

    # handle metadata
    source = "filter_interactions.py"
    datasetMetaProps = {"source": source, "description": "Filter by interactions"}
    clsMappings = {
        # "EnumChargesSrcMolUUID": "java.lang.String",
        # "EnumChargesSrcMolIdx": "java.lang.Integer"
    }
    fieldMetaProps = [
        # {"fieldName": "EnumChargesSrcMolUUID", "values": {"source": source, "description": "UUID of source molecule"}},
        # {"fieldName": "EnumChargesSrcMolIdx", "values": {"source": source, "description": "Index of source molecule"}}
    ]

    input, suppl = rdkit_utils.default_open_input(args.input, args.informat)
    output, writer, output_base = rdkit_utils.default_open_output(args.output,
                                                                 'filter_interactions', args.outformat,
                                                                 thinOutput=False, valueClassMappings=clsMappings,
                                                                 datasetMetaProps=datasetMetaProps,
                                                                 fieldMetaProps=fieldMetaProps,
                                                                 compress=not args.no_gzip)
    report_file = open(output_base + '.report', 'wt')
    count, total, errors = execute(suppl, writer, report_file, args.group_by_field, args.score_field, args.score_descending,
                                   args.stats_fields)

    utils.log(count, total, errors)

    if input:
        input.close()
    writer.flush()
    writer.close()
    output.close()
    report_file.close()

    # re-write the metadata as we now know the size
    if args.outformat == 'json':
        utils.write_squonk_datasetmetadata(output_base, False, clsMappings, datasetMetaProps, fieldMetaProps, size=total)

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': count, '__OutputCount__': total, '__ErrorCount__': errors,
                                          'FilterInteractions': total})


if __name__ == "__main__":
    main()
