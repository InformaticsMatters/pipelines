# Piplelines.

The project experiments with ways to generate data processing piplelines.
The aim is to generate some re-usable building blocks that can be piped 
together into more functional pipelines. Their prime initial use is as executors
for the Squonk Computational Notebook (http://squonk.it) though it is expected
that they might have uses in other environments.

As well as being executable directly they can also be executed in Docker
containers (separately or as a single pipeline). Additionally they can be 
executed using Nextflow (http://nextflow.io) to allow running large jobs 
on HPC-like environments.

Currently it has some python scripts using RDKit (http://rdkit.org) to provide 
basic cheminformatics and comp chem functionality, though other tools will 
be coming soon, including some from the Java ecosystem.

Note: this is highly experimental, everything is subject to change, and 
there are no guarantees that anything works!
That said, if you are interested let me know, and join the fun.

## General principles

### Modularity
Each component should be small but useful. Try to split complex tasks into 
reusable steps. Think how the same steps could be used in other workflows.
Allow parts of ome component to be used in another component where appropriate
but avoid over use. For example see the use of functions in rdkit/conformers.py 
to generate conformers in o3dAlign.py 

### Consistency

Consistent approach to how components function, regarding:

1. Use as simple command line tools that can be piped together
1. Input and outputs either as files of using STDIN and STDOUT
1. Any info/logging written to STDERR to keep STDOUT free for output
1. Consistent approach to command line arguments across components

Generally use consistent coding styles e.g. PEP8 for Python.



## Contact

Any questions contact: 

Tim Dudgeon
tdudgeon@informaticsmatters.com
