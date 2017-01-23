#!/usr/bin/env python

import utils, os
import argparse


### start main execution #########################################

def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit smarts filter')
    utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    args = parser.parse_args()
    utils.log("Screen Args: ", args)

    if not args.output:
        raise ValueError("Must specify output location")

    input, output, suppl, writer, output_base = utils.default_open_input_output(args.input, args.informat, args.output,
                                                                                'smarts_filter', args.outformat)

    ### Define the filter chooser - lots of logic possible
    # SMARTS patterns are defined in poised_filter.py. Currently this is hardcoded.
    # Should make this configurable so that this can be specified by the user at some stage.
    poised_filter = True
    if poised_filter == True:
        from poised_filter import Filter
        filter_to_use = Filter()

    i = 0
    count = 0

    dir_base = os.path.dirname(args.output)

    writer_dict = filter_to_use.get_writers(dir_base)

    for mol in suppl:
        i += 1
        if mol is None: continue
        # Return a dict/class here - indicating which filters passed
        filter_pass = filter_to_use.pass_filter(mol)
        if filter_pass:
            count += 1
            if not args.quiet:
                utils.log(i, "passed")
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            writer.write(mol)
            for reaction in filter_pass:
                writer_dict[reaction].write(mol)
                writer_dict[reaction].flush()

    utils.log("Found", count, " appropriate molecules from "+str(i))
    utils.log("Mols found in: "+ dir_base)

    writer.flush()
    writer.close()
    if input:
        input.close()
    if output:
        output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'SmartsFilter': count})


if __name__ == "__main__":
    main()
