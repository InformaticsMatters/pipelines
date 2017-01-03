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

Tim Dudgeon
tdudgeon@informaticsmatters.com
