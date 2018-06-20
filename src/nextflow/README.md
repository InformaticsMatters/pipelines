# Nextflow Examples

This directory contains examples of running various pipelines components using [Nextflow](http://nextflow.io).

The examples are designed to be run from the top level project directory, and, unless you have Nextflow and all the 
necessary dependencies installed are best run in Docker using the appropriate Docker image. Dockerfiles for these can
be found in the < project root >/src directory and images can be found on Docker Hub.

## File name conventions.

Nextflow files names *.nsd.nf are destined for Squonk (and may have corresponding *.nsd.yml and *.nsd.config files).
For normal use just create files like *.nf. The .nsd.yml extension is needed for it to be recognised by Squonk.

## RDKit implementations

### screen - a simple example

This is a simple one step process that screens an SD file for similarity to a query structure. There's no real need
to use Nextflow for this as it can be run as a single simple command, but its included here as a test that things are 
working.

To run with Docker from the project root dir run this:

```
$ nextflow run src/nextflow/rdkit/screen.nf -c src/nextflow/rdkit/screen.config -with-docker informaticsmatters/rdkit_pipelines
N E X T F L O W  ~  version 0.22.6
Launching `src/nextflow/rdkit/screen.nf` [jovial_hoover] - revision: 1f3a8d73e2
[warm up] executor > local
[3f/31a35a] Submitted process > rdkitScreen (1)
Results: /Users/timbo/dev/git/pipelines/work/3f/31a35a855ec2ab2904f0758aa760a5/results.sdf.gz
```


Those examples use built in defaults for everything. For a more realistic example you would specify the query structure
and the SD file to be screened.

```
nextflow run src/nextflow/rdkit/screen.nf --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --target data/dhfr_3d.sdf.gz --simmin 0.6 -with-docker informaticsmatters/rdkit_pipelines
```    

### screen and conformers - a real pipeline

This is an example that screens an SD file as in the previous example and then generates conformers of the structures 
that passed the screening step. As such it's an example of a two step pipeline.

```
$ nextflow run src/nextflow/rdkit/screen+conformers.nf -c  src/nextflow/rdkit/screen.config -with-docker
N E X T F L O W  ~  version 0.22.6
Launching `src/nextflow/rdkit/screen+conformers.nf` [sad_poincare] - revision: 67201accca
[warm up] executor > local
[f0/191df1] Submitted process > rdkitScreen (1)
[88/72fd4d] Submitted process > rdkitConformer (1)
Results: /Users/timbo/dev/git/pipelines/work/88/72fd4de1f7186f244efd8db7f0d6af/results.sdf.gz
```

Of course the notes in the previous section about running in a Docker container and providing real parameters also apply here.

## PLI scoring

```
$ nextflow run src/nextflow/docking/plip.nf -c src/nextflow/docking/plip.config -with-docker --ligands data/smog/confs.sdf --protein data/smog/DCP2_1.pdb informaticsmatters/pli:latest
N E X T F L O W  ~  version 0.29.0
Launching `src/nextflow/docking/plip.nf` [lethal_northcutt] - revision: 41ef9aebb5
[warm up] executor > local
[6e/6fc729] Submitted process > sdsplit
[77/d0a6a7] Submitted process > pli_scoring (2)
[61/c89746] Submitted process > pli_scoring (1)
[2e/39b06e] Submitted process > pli_scoring (4)
[6b/e59ad2] Submitted process > pli_scoring (3)
[95/6ea8c9] Submitted process > pli_scoring (9)
[8a/16f929] Submitted process > pli_scoring (8)
[c8/0ebfcf] Submitted process > pli_scoring (10)
[d4/4da278] Submitted process > pli_scoring (6)
[0a/fda450] Submitted process > pli_scoring (7)
[60/09e0d0] Submitted process > pli_scoring (11)
[19/1c7ca9] Submitted process > pli_scoring (5)
[58/720b16] Submitted process > pli_scoring (12)
[4e/2524ad] Submitted process > pli_scoring (13)
[dc/11f9b9] Submitted process > pli_scoring (14)
[98/8a0523] Submitted process > pli_scoring (15)
[a6/0b1ad3] Submitted process > results
```

## SMoG2016 scoring

```
$ nextflow run src/nextflow/docking/smog.nf -c src/nextflow/docking/smog.config -with-docker --ligands data/smog/confs.sdf --protein data/smog/DCP2_1.pdb informaticsmatters/smog:latest
N E X T F L O W  ~  version 0.29.0
Launching `src/nextflow/docking/smog.nf` [cheesy_dijkstra] - revision: 220c72f4e2
[warm up] executor > local
[70/2ca9ff] Submitted process > sdsplit
[a9/b85a58] Submitted process > smog_scoring (3)
[56/1e0219] Submitted process > smog_scoring (2)
[82/749fa3] Submitted process > smog_scoring (4)
[b5/8062b0] Submitted process > smog_scoring (1)
[a5/225732] Submitted process > smog_scoring (8)
[76/c34cd8] Submitted process > smog_scoring (10)
[0e/4e5449] Submitted process > smog_scoring (7)
[90/e905c3] Submitted process > smog_scoring (12)
[db/c980d1] Submitted process > smog_scoring (6)
[77/0cfb16] Submitted process > smog_scoring (9)
[77/5fcc28] Submitted process > smog_scoring (5)
[ee/6a7f68] Submitted process > smog_scoring (11)
[77/5b087f] Submitted process > smog_scoring (13)
[61/9b0853] Submitted process > smog_scoring (15)
[70/5efdec] Submitted process > smog_scoring (14)
[a6/06febe] Submitted process > results
```

## rDock docking

```
$ nextflow run src/nextflow/docking/rdock.nf -c src/nextflow/docking/rdock.config -with-docker --ligands data/hivpr_ligprep_100.sdf.gz --protein data/hivpr_rdock.mol2 --asfile data/hivpr_rdock.as --prmfile data/hivpr_rdock.prm --num_dockings 10
N E X T F L O W  ~  version 0.29.0
Launching `src/nextflow/docking/rdock.nf` [festering_newton] - revision: e2ebfb5346
[warm up] executor > local
[51/31d685] Submitted process > sdsplit
[db/e784c2] Submitted process > rdock (3)
[aa/9927f7] Submitted process > rdock (1)
[74/7b5c3b] Submitted process > rdock (2)
[e1/c3301c] Submitted process > rdock (4)
[bc/87e17d] Submitted process > results
$
```


## Usage in Squonk

To be used in Squonk you must produce a `*.nsd.yml` file with the Squonk Service Descriptor. This has a corresponding Nextflow
file and configuration. These can either be specified directly in the `*.nsd.yml` file, or as separate files named `*.nsd.nf` and 
`*.nsd.config` (using the same base file name as the .yml file).
See [src/nextflow/screen-dataset.nsd.yml] for an example.

When Squonk executes your pipeline it does the following:
* creates a temporary working directory
* copies your input files (names specified in the service descriptor) to this directory
* copies your `*.nsd.nf` file to the file `nextflow.nf` in this directory
* copies your `*.nsd.conf` file to the file `nextflow.conf` in this directory
* creates a simple `execute` bash script that will executes the pipeline with the options provided by the user and as specified 
in the service descriptor
* Runs a container using the `informaticsmatters/nextflow` Docker image, mounting in the working directory and running the
`execute` script
* After execution grabs the expected output files and as the results of execution
* deletes the temporary working directory (unless working in debug mode)

To test basic execution as it would take place in Squonk create a directory containing the input files that are needed by
the pipeline (make sure they are named as expected by the pipeline) and then execute using the `execute` shell script that 
would be generated by squonk. This file has a very simple content, something like this:
```
#!/bin/sh
nextflow run nextflow.nf -c nextflow.config <your-parameters> -with-docker
```
To explain this a bit:
* `nextflow.nf` is your `*.nsd.nf` file
* `-c nextflow.config` specifies which config file to use (your `*.nsd.nf` file) 
* `<your-parameters>` is what is generated from the `nextflowParams` property in your `*.nsd.yaml` file.

An example of how to execute would be this:

```
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow sh -c 'nextflow run plip.nsd.nf -c plip.nsd.config --score 100.0 -with-docker'
``` 

See [test-nextflow.sh](../..test-nextflow.sh) for some real examples.
