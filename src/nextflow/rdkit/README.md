# Nextflow Examples

This directory contains examples of running various pipelines components using [Nextflow](http://nextflow.io).

The examples are designed to be run from the top level project directory, and, unless you have Nextflow and all the 
necessary dependencies installed are best run in Docker using the appropriate Docker image. Dockerfiles for these can
be found in the < project root >/src directory and images can be found on Docker Hub.

## File name conventions.

Nextflow files names *.nsd.nf are destined for Squonk (and may have corresponding *.nsd.yml and *.nsd.config files).
For normal use just create files like *.nf. The .nsd.yml extension is needed for it to be recognised by Squonk.

## RDKit examples

### screen - a simple example

This is a simple one step process that screens an SD file for similarity to a query structure. There's no real need
to use Nextflow for this as it can be run as a single simple command, but its included here as a test that things are 
working.

To run with Docker from the project root dir run this:

```sh
$ nextflow run src/nextflow/rdkit/screen.nf -with-docker informaticsmatters/rdkit_pipelines
N E X T F L O W  ~  version 0.22.6
Launching `src/nextflow/rdkit/screen.nf` [jovial_hoover] - revision: 1f3a8d73e2
[warm up] executor > local
[3f/31a35a] Submitted process > rdkitScreen (1)
Results: /Users/timbo/dev/git/pipelines/work/3f/31a35a855ec2ab2904f0758aa760a5/results.sdf.gz
```


Those examples use built in defaults for everything. For a more realistic example you would specify the query structure
and the SD file to be screened.

```sh
nextflow run src/nextflow/rdkit/screen.nf --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --target data/dhfr_3d.sdf.gz --simmin 0.6 -with-docker informaticsmatters/rdkit_pipelines
```    

### screen and conformers - a real pipeline

This is an example that screens an SD file as in the previous example and then generates conformers of the structures 
that passed the screening step. As such it's an example of a two step pipeline.

```sh
$ nextflow run src/nextflow/rdkit/screen+conformers.nf -with-docker informaticsmatters/rdkit_pipelines
N E X T F L O W  ~  version 0.22.6
Launching `src/nextflow/rdkit/screen+conformers.nf` [sad_poincare] - revision: 67201accca
[warm up] executor > local
[f0/191df1] Submitted process > rdkitScreen (1)
[88/72fd4d] Submitted process > rdkitConformer (1)
Results: /Users/timbo/dev/git/pipelines/work/88/72fd4de1f7186f244efd8db7f0d6af/results.sdf.gz
```

Of course the notes in the previous section about running in a Docker container and providing real parameters also apply here.

## Usage in Squonk

To be used in Squonk you must produce a *.nsd.yml file with the Squonk Service Descriptor. This has a corresponding Nextlow
file and configuration. These can either be specified directly in the .yml file, or as separate files named *.nsd.nf and 
*.nsd.config (using the same base file name as the .yml file).
See [src/nextflow/screen-dataset.nsd.yml] for an example.