#!/usr/bin/env python

"""The "Squonk Cell Design Template"

Use this template and fill out the design sections to submit a cell
implementation along with a detailed description of
the cell's behaviour for design.

With regard to the implementation, the code, you can supply your that
separately if that's easier for you or you can insert your implementation
into the 'implementation()' function. If you insert code into this template
don't worry about syntax errors or successful compilation at this
stage. The template is simply a medium you need to use to submit a request
for a new cell implementation that we can implement and test for you.

By filling out the sections you'll find below you'll helps us to: -

-   Understand the detailed requirements of the cell
-   The cell's parameters and its input and output files
-   Test the cell using parameter values and/or test files that you provide

So that we can build and test your cell please fill out the following design
sections: -

#-----------#
# Section A # -- IDENTITY -----------------------------------------------------
#-----------#

Use this section to provide a unique identity for your cell.
This is used in Squonk for cell identification, metrics collection and
searching.

You need to provide the following fields that wil be used to
distinguish your cell in Squonk: -

    - name: <compact, simple phrase description>
      description: <description>
      python: <minimum..maximum python versions>
      tags:
      - <search term>
      icon: <name of PNG icon file, supplied separately|-na->

If a value is not applicable use the text '-na-' for the parameter value.

The following example should help to illustrate how these fields
might be defined: -

- name: NOAEL Calculator
  description: A NOAEL (No Observed Adverse Effect Level) calculator
  python: 3.8..
  tags:
  - noael
  icon: -na-

Insert your identity definition, just like the example above,
between the following IDENTITY BEGIN and IDENTITY END markers: -

: IDENTITY BEGIN
: IDENTITY END

#-----------#
# Section B # -- INPUTS (Parameters) ------------------------------------------
#-----------#

Use this section to provide a list of descriptions for any user-defined
parameters that control your cell's behaviour. These are typically the
command-line parameters or code variables that that you use for your
implementation and you expect a user to legitimately modify. You should
list every parameter or variable here. This list should not
include internal variables or constants, only parameters or variables
you expect a user to change.

You need to provide the following fields for each variable you want
to expose as an input to the cell: -

    - label: <compact symbolic name>
      description: <compact descriptive phrase>
      type: <string|integer|float>
      required: <yes|no>
      range: <min..max|choice>
      default_value: <value|-na->
      dependency: <dependency on any other variable's value|-na->

The following example should help to illustrate how these fields
might be defined. here we have a numerical input parameter that must be
provided by the user (i.e. has no default value). It has a minimum value of
1 but no upper limit on the value, has no default value and is not
dependent on any other user-defined cell variable: -

- label: NOAEL
  description: No Observable Adverse Effect Level (mg/kg/day)
  type: integer
  required: yes
  range: 1..
  default_value: -na-
  dependency: -na-

If you have any variables then insert a parameter definition for each one,
just like the example above, between the following PARAMETER DEFINITIONS BEGIN
and PARAMETER DEFINITIONS END markers: -

: PARAMETER DEFINITIONS BEGIN
: PARAMETER DEFINITIONS END

#-----------#
# Section C # -- INPUTS (Files) -----------------------------------------------
#-----------#

If your cell's behaviour is driven by files what types of file to you expect
to be able to attach to the cell and should we be aware of an expectation to
handle very large files (e.g. greater than, say, 50MB)?

For each file provide the following: -

    - type: <type>
      special_notes: <any notes or special fields that are important>
      maximum_size: <size>

The following example tells us that you expect to be presented with
an SDF file (compressed or uncompressed) no larger than 5MB: -


- type: sdf (compressed or uncompressed)
  special_notes: -na-
  maximum_size: 5M

If you expect to process any files then insert a file definition for each one,
just like the example above, between the following INPUT FILES BEGIN
and INPUT FILES END markers: -

: INPUT FILES BEGIN
: INPUT FILES END

#-----------#
# Section D # -- OUTPUTS ------------------------------------------------------
#-----------#

Here we need to know what you want to 'emit' from your cell, i.e. what does
it produce? This could be as simple as a text message like "Calculation is 5.6"
or something more complicated like a table of results, a file
(like an SDF file) or an image (a matplot-lib-like plot).

Provide at least one output definition for each of your cell's outputs,
just like the example above, between the following OUTPUTS BEGIN and
OUTPUTS END markers: -

: OUTPUTS BEGIN
: OUTPUTS END

#-----------#
# Section E # -- IMPLEMENTATION (REQUIREMENTS) --------------------------------
#-----------#

If your code has dependencies (requirements) on external Python modules
that can be found on PyPI, for example 'numpy' or 'pillow', then please
identify them here.

We strongly encourage you to be very explicit with regard to your module
versions, especially by providing a minimum version. Just relying
on a module by name does not protect your implementation from breakages
if a new (major) version of the module is published.

Paste your requirements, typically the content of your requirements file
between the following REQUIREMENTS BEGIN and REQUIREMENTS END markers: -

: REQUIREMENTS BEGIN
: REQUIREMENTS END

#-----------#
# Section F # -- CELL TESTS ---------------------------------------------------
#-----------#

In the preceding sections you've defined your input parameters, input files
output files and module dependencies. Now, to ensure your cell behaves as
expected you need to tell us what you expect the outputs to look like when
the cell is run.

To do this you must provide us with a set of tests using example parameter
values and input files (where appropriate).

You must provide at least one test.

If you have test data as well as parameter values you must also let us have
those.

For each test provide the following: -

    - name: <compact identifying name for the test>
      parameters:
      - <name>: <value>
      files:
      - <name>
      status_message: <description|-na->
      output_file: <description|-na->

The following example tells us that the test is based on the input file
'example.sdf' using '56' as the 'NOAEL' parameter value. The result is
a message in the cell status line with the text "Result is 4.5" and there is no
output file: -

- name: Test NOEL 56
  parameters:
  - NOAEL: 56
  files:
  - example.sdf
  status_message: Result is 4.5
  output_file: -na-

Insert test definitions, just like the example above, between the following
CELL TESTS BEGIN and CELL TESTS END markers: -

: CELL TESTS BEGIN
: CELL TESTS END

#-----------#
# Section G # -- ANYTHING ELSE ------------------------------------------------
#-----------#

Finally, is there anything else you think we need to know about the
implementation or the behaviour of your cell that has not already been covered
by the preceding sections?

Provide extra notes between the following EXTRA BEGIN and EXTRA END markers: -

: EXTRA BEGIN
: EXTRA END

"""
import argparse
import os

from pipelines_utils import parameter_utils, utils

_CELL_NAME: str = 'Template'
_CELL_KEY: str = 'template'


def implementation() -> None:
    # CUT FROM HERE...
    #
    # REPLACE ME AND THE pass LINE WITH YOUR CODE FRAGMENT
    # DO NOT WORRY ABOUT SYNTAX OR EDITOR WARNINGS AND ERRORS
    # RELATING TO YOUR CODE. THE CODE IS NOT EXPECTED TO WORK
    # SIMPLY BY PASTING IT HERE - THE CODE, COMBINED WITH THE
    # ANSWERS TO THE QUESTIONS IN THE MODULE'S DOCUMENTATION BLOCK
    # WILL HELP US UNDERSTAND HOW TO CONNECT IT TO SQUONK'S
    # EXECUTION ENVIRONMENT.
    pass
    # ...TO HERE


def main() -> None:

    # Command line argument definitions #######################################

    parser = argparse.ArgumentParser(description=_CELL_NAME)
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    args = parser.parse_args()
    utils.log("Args: ", args)


if __name__ == "__main__":
    main()
