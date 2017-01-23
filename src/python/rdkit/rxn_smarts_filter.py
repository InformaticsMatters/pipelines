#!/usr/bin/env python

import utils, os
import argparse


### start main execution #########################################

def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit rxn smarts filter')
    utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-m', '--multi', action='store_true', help='Output one file for each reaction')

    args = parser.parse_args()
    utils.log("Screen Args: ", args)

    if not args.output and args.multi:
        raise ValueError("Must specify output location when writing individual result files")

    input, output, suppl, writer, output_base = utils.default_open_input_output(args.input, args.informat, args.output,
                                                                                'rxn_smarts_filter', args.outformat)

    ### Define the filter chooser - lots of logic possible
    # SMARTS patterns are defined in poised_filter.py. Currently this is hardcoded.
    # Should make this configurable so that this can be specified by the user at some stage.
    poised_filter = True
    if poised_filter == True:
        from poised_filter import Filter
        filter_to_use = Filter()

    i = 0
    count = 0



    if args.multi:
        dir_base = os.path.dirname(args.output)
        writer_dict = filter_to_use.get_writers(dir_base)
    else:
        writer_dict = None
        dir_base = None

    for mol in suppl:
        i += 1
        if mol is None: continue
        # Return a dict/class here - indicating which filters passed
        filter_pass = filter_to_use.pass_filter(mol)
        utils.log("Found",str(len(filter_pass)),"matches")

        if filter_pass:
            count += 1

            for reaction in filter_pass:
                # TODO - write something more useful here.
                # Ideally the part of the molecule that matched so that it can be highlighted
                # But note that there can potentially be multiple matches
                mol.SetProp(reaction, "pass")
                if args.multi:
                    writer_dict[reaction].write(mol)
                    writer_dict[reaction].flush()

            writer.write(mol)
            writer.flush()


    utils.log("Matched", count, "molecules from a total of", i)
    if dir_base:
        utils.log("Individual SD files found in: "+ dir_base)

    writer.flush()
    writer.close()
    if input:
        input.close()
    if output:
        output.close()
    # close the individual writers
    if writer_dict:
        for key in writer_dict:
            writer_dict[key].close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RxnSmartsFilter': count})


if __name__ == "__main__":
    main()
