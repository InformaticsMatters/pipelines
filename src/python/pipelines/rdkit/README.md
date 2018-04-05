# RDKit examples

See test-python.sh for working examples. The ones listed here could be out of date.

Examples of running some of these examples with the [Nextflow](http://nextflow.io) workflow system
can be found [here](../../../nextflow/rdkit/README.md)

Run these commands for the root directory of this project.

Install the python modules using (from the top levlel directory of this project):

```sh
pip install -e src/python
Obtaining file:///Users/timbo/dev/git/pipelines/src/python
Installing collected packages: pipelines
  Found existing installation: pipelines 0.0.1b0
    Uninstalling pipelines-0.0.1b0:
      Successfully uninstalled pipelines-0.0.1b0
  Running setup.py develop for pipelines
Successfully installed pipelines
```


## 1. screen.py
Performs ligand based virtual screening using RDKit descriptors

### Native execution

These examples need Python and RDKit to be present

```sh
python -m pipelines.rdkit.screen -h
usage: screen.py [-h] [--qsmiles QSMILES | --qmolfile QMOLFILE]
                 [--simmin SIMMIN] [--simmax SIMMAX]
                 [-d {rdkit,morgan3,morgan2,maccs}]
                 [-m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}]
                 [-f {hac,mw}] [--hacmin HACMIN] [--hacmax HACMAX]
                 [--mwmin MWMIN] [--mwmax MWMAX] [-i INPUT] [-if {sdf,json}]
                 [-o OUTPUT] [-of {sdf,json}] [--meta] [--thin] [-q]

RDKit screen

optional arguments:
  -h, --help            show this help message and exit
  --qsmiles QSMILES     query structure as smiles (incompatible with -qmolfile
                        arg)
  --qmolfile QMOLFILE   query structure as filename in molfile format
                        (incompatible with -qsmiles arg)
  --simmin SIMMIN       similarity lower cutoff (1.0 means identical)
  --simmax SIMMAX       similarity upper cutoff (1.0 means identical)
  -d {rdkit,morgan3,morgan2,maccs}, --descriptor {rdkit,morgan3,morgan2,maccs}
                        descriptor or fingerprint type (default rdkit)
  -m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}, --metric {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}
                        similarity metric (default tanimoto)
  -f {hac,mw}, --fragment {hac,mw}
                        Find single fragment if more than one (hac = biggest
                        by heavy atom count, mw = biggest by mol weight )
  --hacmin HACMIN       Min heavy atom count
  --hacmax HACMAX       Max heavy atom count
  --mwmin MWMIN         Min mol weight
  --mwmax MWMAX         Max mol weight
  -i INPUT, --input INPUT
                        Input file, if not defined the STDIN is used
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
  --thin                Thin output mode
  -q, --quiet           Quiet mode
```

```sh
gunzip -c data/dhfr_3d.sdf.gz | python -m pipelines.rdkit.screen --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --simmin 0.5 --informat sdf | gzip > tmp/screened_mols.sdf.gz
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed smiles. Gzip the resulting SDF and write to screened_mols.sdf.gz

```sh
python -m pipelines.rdkit.screen --qmolfile data/pyrimethamine.mol -i data/dhfr_3d.sdf.gz --simmin 0.5 --simmax 0.8 -if sdf -of sdf > screened_mols.sdf
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed file in molfile format of between 0.5 and 0.8. Write results to screened_mols.sdf

### Execution in Docker

From the rdkit directory (adjust arguments if different):

```sh
gunzip -c data/dhfr_3d.sdf.gz | docker run -i --rm -v $PWD:/work -w /work informaticsmatters/rdkit_pipelines python -m pipelines.rdkit.screen --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --simmin 0.49 -if sdf -o tmp/myfile.sdf
Screen Args:  Namespace(descriptor='rdkit', fragment=None, hacmax=None, hacmin=None, informat='sdf', input=None, meta=False, metric='tanimoto', mwmax=None, mwmin=None, outformat=None, output='tmp/myfile.sdf', qmolfile=None, qsmiles='C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2', quiet=False, simmax=1.0, simmin=0.49, thin=False)
No output format specified - using sdf
718 0.49537037037
723 0.490173410405
745 0.566194837635
Found 3 similar molecules
```
Input pipe into Docker container using STDIN.
Results sent from Docker container to STDOUT and redirected to file tmp/myfile.sdf

The Docker image informaticsmatters/rdkit_pipelines can be found 
[here](https://hub.docker.com/r/informaticsmatters/rdkit_pipelines/) 
and is automatically updated whenever this Git repository is updated.


## 2. conformers.py
Generates 3D conformers

```sh
python -m pipelines.rdkit.conformers -h
usage: conformers.py [-h] [-n NUM] [-a ATTEMPTS] [-r RMSD] [-c {rmsd,tdf}]
                     [-t THRESHOLD] [-e EMIN] [-i INPUT] [-if {sdf,json}]
                     [-o OUTPUT] [-of {sdf,json}] [--meta] [--smiles SMILES]

RDKit conformers

optional arguments:
  -h, --help            show this help message and exit
  -n NUM, --num NUM     number of conformers to generate
  -a ATTEMPTS, --attempts ATTEMPTS
                        number of attempts
  -r RMSD, --rmsd RMSD  prune RMSD threshold
  -c {rmsd,tdf}, --cluster {rmsd,tdf}
                        Cluster method (RMSD or TFD). If None then no
                        clustering
  -t THRESHOLD, --threshold THRESHOLD
                        cluster threshold (default of 2.0 for RMSD and 0.3 for
                        TFD)
  -e EMIN, --emin EMIN  energy minimisation iterations (default of 0 means
                        none)
  -i INPUT, --input INPUT
                        Input file, if not defined the STDIN is used
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
  --smiles SMILES       input structure as smiles (incompatible with using
                        files or stdin for input)
```

```sh
gunzip -c data/Kinase_inhibs.sdf.gz | python -m pipelines.rdkit.conformers -n 5 -if sdf > out.sdf
```
Tries to generate 5 conformers for each structure in the input file (may generate less if 2 conformers differ by less than the default RMSD) 

```sh
python -m pipelines.rdkit.conformers -i data/Kinase_inhibs.sdf.gz -n 5 -a 20 -c rmsd > out.sdf
```
Input specified using -i argument. Tries to generate 5 conformers using a maximum of 20 attempts, and clusters the conformers using RMSD distance.


## 3. o3dAlign.py
Performs rigid 3D alignment of structures
```sh
python -m pipelines.rdkit.o3dAlign -h
usage: o3dAlign.py [-h] [--qmolidx QMOLIDX] [-t THRESHOLD] [-n NUM]
                   [-a ATTEMPTS] [-r RMSD] [-e EMIN] [-i INPUT]
                   [-if {sdf,json}] [-o OUTPUT] [-of {sdf,json}] [--meta]
                   query

Open3DAlign with RDKit

positional arguments:
  query                 query molfile

optional arguments:
  -h, --help            show this help message and exit
  --qmolidx QMOLIDX     Query molecule index in SD file if not the first
  -t THRESHOLD, --threshold THRESHOLD
                        score cuttoff relative to alignment of query to itself
  -n NUM, --num NUM     number of conformers to generate, if None then input
                        structures are assumed to already be 3D
  -a ATTEMPTS, --attempts ATTEMPTS
                        number of attempts to generate conformers
  -r RMSD, --rmsd RMSD  prune RMSD threshold for excluding conformers
  -e EMIN, --emin EMIN  energy minimisation iterations for generated confomers
                        (default of 0 means none)
  -i INPUT, --input INPUT
                        Input file, if not defined the STDIN is used
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
```


```sh
gunzip -c data/dhfr_3d.sdf.gz | python -m pipelines.rdkit.o3dAlign data/pyrimethamine.mol -n 5 -t 10 -if sdf > my_output.sdf
```
Aligns structures from STDIN to the specified molecule (3D), generating up to 5 conformers to align. Redirects output structures to my_output.sdf.

```sh
gunzip -c data/dhfr_3d.sdf.gz | python -m pipelines.rdkit.o3dAlign data/pyrimethamine.mol -if sdf -n 5 -t 10 -o screen1
```

```sh
python o3dAlign.py ../data/pyrimethamine.mol -i ../data/dhfr_3d.sdf.gz -o o3da
```
Aligns the 3D structures in the file dhfr_3d.sdf.gz to the 3D structure pyrimethamine.mol writing the output to o3da.sdf.gz

## 4. filter.py
Filters and fixes molecules
```sh
python -m pipelines_utils_rdkit.filter -h
usage: filter.py [-h] [-f {hac,mw}] [--hacmin HACMIN] [--hacmax HACMAX]
                 [--mwmin MWMIN] [--mwmax MWMAX] [-l LIMIT] [-c CHUNKSIZE]
                 [-d DIGITS] [--no-gzip] [--thin] [-q] [-i INPUT]
                 [-if {sdf,json}] [-o OUTPUT] [-of {sdf,json}] [--meta]

RDKit filter

optional arguments:
  -h, --help            show this help message and exit
  -f {hac,mw}, --fragment {hac,mw}
                        Find single fragment if more than one (hac = biggest
                        by heavy atom count, mw = biggest by mol weight )
  --hacmin HACMIN       Min heavy atom count
  --hacmax HACMAX       Max heavy atom count
  --mwmin MWMIN         Min mol weight
  --mwmax MWMAX         Max mol weight
  -l LIMIT, --limit LIMIT
                        Limit output to this many records
  -c CHUNKSIZE, --chunksize CHUNKSIZE
                        Split output into chunks of size c. Output will always
                        be files. Names like filter1.sdf.gz, filter2.sdf.gz
                        ...
  -d DIGITS, --digits DIGITS
                        When splitting zero pad the file name to this many
                        digits so that they are in sorted order. Names like
                        filter001.sdf.gz, filter002.sdf.gz ...
  --no-gzip             Do not compress the output (STDOUT is never compressed
  --thin                Thin output mode
  -q, --quiet           Quiet mode - suppress reporting reason for filtering
  -i INPUT, --input INPUT
                        Input file, if not defined the STDIN is used
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
```

```sh
gunzip -c data/Kinase_inhibs.sdf.gz | python -m pipelines_utils_rdkit.filter --hacmin 25 --hacmax 30 -if sdf
```


## 5. cluster_butina.py
Molecule clustering using Butina method
```sh
python -m pipelines.rdkit.cluster_butina -h
usage: cluster_butina.py [-h] [-t THRESHOLD]
                         [-d {rdkit,morgan3,morgan2,maccs}]
                         [-m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}]
                         [-q] [-n NUM] [-e EXCLUDE] [-f FIELD] [--min | --max]
                         [-i INPUT] [-if {sdf,json}] [-o OUTPUT]
                         [-of {sdf,json}] [--meta] [--thin]

RDKit Butina Cluster

optional arguments:
  -h, --help            show this help message and exit
  -t THRESHOLD, --threshold THRESHOLD
                        similarity clustering threshold (1.0 means identical)
  -d {rdkit,morgan3,morgan2,maccs}, --descriptor {rdkit,morgan3,morgan2,maccs}
                        descriptor or fingerprint type (default rdkit)
  -m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}, --metric {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}
                        similarity metric (default tanimoto)
  -q, --quiet           Quiet mode
  -n NUM, --num NUM     maximum number to pick for diverse subset selection
  -e EXCLUDE, --exclude EXCLUDE
                        threshold for excluding structures in diverse subset
                        selection (1.0 means identical)
  -f FIELD, --field FIELD
                        field to use to optimise diverse subset selection
  --min                 pick lowest value specified by the --field option
  --max                 pick highest value specified by the --field option
  -i INPUT, --input INPUT
                        Input file, if not defined the STDIN is used
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
  --thin                Thin output mode
```

```sh
gunzip -c data/Kinase_inhibs.sdf.gz | python -m pipelines.rdkit.cluster_butina -t 0.6 -if sdf
```

## 6. Connecting modules using pipes

Pipelines modules are designed to be connected using pipes.
 
```sh
python -m pipelines.rdkit.screen --qmolfile data/pyrimethamine.mol -i data/dhfr_3d.sdf.gz --simmin 0.5 --simmax 0.8 -if sdf -of sdf |\
 python -m pipelines.rdkit.conformers -n 5 -if sdf |\
 gzip > tmp/my_output.sdf.gz
```
Screens a set of compounds for similarity, generates conformers for those that pass, 
gzip the results and write to file.
