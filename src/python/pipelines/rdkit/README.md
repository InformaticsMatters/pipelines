#RDKit examples

See test-python.sh for working examples. The ones listed here could be out of date.

##1. screen.py
Performs ligand based virtual screening using RDKit descriptors

### Native execution

These examples need Python and RDKit to be present

```sh
usage: screen.py [-h] [--smiles SMILES] [--molfile MOLFILE] [--simmin SIMMIN]
                 [--simmax SIMMAX] [-d {maccs,morgan2,morgan3,rdkit}]
                 [-m {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}]
                 [-f {hac,mw}] [--hacmin HACMIN] [--hacmax HACMAX]
                 [--mwmin MWMIN] [--mwmax MWMAX] [-i INPUT] [-o OUTPUT]
                 [-if {sdf,json}] [-of {sdf,json}] [-q]

RDKit screen

optional arguments:
  -h, --help            show this help message and exit
  --smiles SMILES       query structure as smiles (incompatible with -molfile
                        arg)
  --molfile MOLFILE     query structure as filename in molfile format
                        (incompatible with -smiles arg)
  --simmin SIMMIN       similarity lower cutoff (1.0 means identical)
  --simmax SIMMAX       similarity upper cutoff (1.0 means identical)
  -d {maccs,morgan2,morgan3,rdkit}, --descriptor {maccs,morgan2,morgan3,rdkit}
                        descriptor or fingerprint type (default rdkit)
  -m {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}, --metric {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}
                        similarity metric (default tanimoto)
  -f {hac,mw}, --fragment {hac,mw}
                        Find single fragment if more than one (hac = biggest
                        by heavy atom count, mw = biggest by mol weight )
  --hacmin HACMIN       Min heavy atom count
  --hacmax HACMAX       Max heavy atom count
  --mwmin MWMIN         Min mol weight
  --mwmax MWMAX         Max mol weight
  -i INPUT, --input INPUT
                        Input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  -q, --quiet           Quiet mode
```

```sh
zcat ../data/dhfr_3d.sdf.gz | python screen.py --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --simmin 0.5 --informat sdf | gzip > screened_mols.sdf.gz
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed smiles. Gzip the resulting SDF and write to screened_mols.sdf.gz

```sh
zcat ../data/dhfr_3d.sdf.gz | python screen.py --molfile ../data/pyrimethamine.mol --simmin 0.5 --simmax 0.8 --informat sdf > screened_mols.sdf
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed file in molfile format of between 0.5 and 0.8. Write results to screened_mols.sdf

### Execution in Docker

From the rdkit directory (adjust arguments if different):

```sh
gunzip -c ../data/dhfr_3d.sdf.gz | docker run -i --rm -v $(pwd):/work -w /work informaticsmatters/rdkit python screen.py --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' --simmin 0.49 --informat sdf > myfile.sdf
Screen Args:  Namespace(descriptor='rdkit', input=None, metric='tanimoto', molfile=None, output=None, simmax=999, simmin=0.49, smiles='C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2')
718 0.49537037037
723 0.490173410405
745 0.566194837635
Found 3 similar molecules
```
Input pipe into Docker container using STDIN.
Results sent from Docker container to STDOUT and redirected to file myfile.sdf

### Execution in Nextflow
Requires Nextflow to be installed.
#### Native
Requires RDKit to be installed


#### Docker
Uses RDKit provided by the container

```sh
nextflow run screen.nf --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' -with-docker informaticsmatters/rdkit
N E X T F L O W  ~  version 0.22.6
Launching `screen.nf` [desperate_saha] - revision: b2fbd583dc
[warm up] executor > local
[88/18d787] Submitted process > rdkitScreen (1)
Results: /Users/timbo/dev/git/pipelines/rdkit/work/88/18d787eea01ff4fd1c15db9e60ff4e/results.sdf.gz
```

## 2. conformers.py
Generates 3D conformers

```sh
python conformers.py -h
usage: conformers.py [-h] [-n NUM] [-a ATTEMPTS] [-r RMSD] [-c {rmsd,tfd}]
                     [-t THRESHOLD] [-e EMIN] [-i INPUT] [-o OUTPUT]
                     [-if {sdf,json}]

RDKit conformers

optional arguments:
  -h, --help            show this help message and exit
  -n NUM, --num NUM     number of conformers to generate
  -a ATTEMPTS, --attempts ATTEMPTS
                        number of attempts
  -r RMSD, --rmsd RMSD  prune RMSD threshold
  -c {rmsd,tfd}, --cluster {rmsd,tfd}
                        Cluster method (rmsd or tfd). If None then no
                        clustering
  -t THRESHOLD, --threshold THRESHOLD
                        cluster threshold (default of 2.0 for RMSD and 0.3 for
                        TFD)
  -e EMIN, --emin EMIN  energy minimisation iterations (default of 0 means
                        none)
  -i INPUT, --input INPUT
                        Input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN as input this needs to
                        be specified.
```

```sh
zcat ../data/Kinase_inhibs.sdf.gz | python conformers.py -n 5 -if sdf > out.sdf
```
Tries to generate 5 conformers for each structure in the input file (may generate less if 2 conformers differ by less than the default RMSD) 

```sh
python conformers.py -i ../data/Kinase_inhibs.sdf.gz -n 5 -a 20 -c rmsd > out.sdf
```
Input specified using -i argument. Tries to generate 5 conformers using a maximum of 20 attempts, and clusters the conformers using RMSD distance.


## 3. o3dAlign.py
Performs rigid 3D alignment of structures
```sh
python o3dAlign.py -h
usage: o3dAlign.py [-h] [-t THRESHOLD] [-n NUM] [-a ATTEMPTS] [-r RMSD]
                   [-e EMIN] [-i INPUT] [-o OUTPUT] [-if {sdf,json}]
                   query

Open3DAlign with RDKit

positional arguments:
  query                 query molfile

optional arguments:
  -h, --help            show this help message and exit
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
                        Input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN as input this needs to
                        be specified.

```


```sh
zcat ../data/dhfr_3d.sdf.gz | python o3dAlign.py ../data/pyrimethamine.mol -n 5 -t 10 -if sdf > my_output.sdf
```
Aligns structures from STDIN to the specified molecule (3D), generating up to 5 conformers to align. Redirects output structures to my_output.sdf.

```sh
zcat ../data/dhfr_3d.sdf.gz | python o3dAlign.py ../data/pyrimethamine.mol -n 5 -t 10 -o screen1
```

```sh
python o3dAlign.py ../data/pyrimethamine.mol -i ../data/dhfr_3d.sdf.gz -o o3da
```
Aligns the 3D structures in the file dhfr_3d.sdf.gz to the 3D structure pyrimethamine.mol writing the output to o3da.sdf.gz

## 4. filter.py
Filters and fixes molecules
```sh
python filter.py -h
usage: filter.py [-h] [-f {hac,mw}] [--hacmin HACMIN] [--hacmax HACMAX]
                 [--mwmin MWMIN] [--mwmax MWMAX] [-l LIMIT] [-c CHUNKSIZE]
                 [-q] [-i INPUT] [-o OUTPUT] [-if {sdf,json}]

RDKit filter

optional arguments:
  -h, --help            show this help message and exit
  -f {hac,mw}, --fragment {hac,mw}
                        Find single fragment if more than one (hac = biggest
                        by heavy atome count, mw = biggest by mol weight )
  --hacmin HACMIN       Min heavy atom count
  --hacmax HACMAX       Max heavy atom count
  --mwmin MWMIN         Min mol weight
  --mwmax MWMAX         Max mol weight
  -l LIMIT, --limit LIMIT
                        Limit output to this many records
  -c CHUNKSIZE, --chunksize CHUNKSIZE
                        Split output into chunks of size c. Output will always
                        be files. Names like filter01.sdf.gz, filter02.sdf.gz
                        ...
  -q, --quiet           Quiet mode - suppress reporting reason for filtering
  -i INPUT, --input INPUT
                        Input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
```

```sh
zcat ../data/Kinase_inhibs.sdf.gz | python filter.py --hacmin 25 --hacmax 30 -if sdf
```


## 5. cluster_butina.py
Molecule clustering using Butina method
```sh
python cluster_butina.py -h
usage: cluster_butina.py [-h] [-t THRESHOLD]
                         [-d {rdkit,morgan3,morgan2,maccs}]
                         [-m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}]
                         [-i INPUT] [-o OUTPUT] [-if {sdf,json}]
                         [-of {sdf,json}] [--meta]

RDKit screen

optional arguments:
  -h, --help            show this help message and exit
  -t THRESHOLD, --threshold THRESHOLD
                        similarity clustering threshold (1.0 means identical)
  -d {rdkit,morgan3,morgan2,maccs}, --descriptor {rdkit,morgan3,morgan2,maccs}
                        descriptor or fingerprint type (default rdkit)
  -m {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}, --metric {sokal,cosine,rogotgoldberg,dice,braunblanquet,asymmetric,kulczynski,mcconnaughey,russel,tanimoto}
                        similarity metric (default tanimoto)
  -i INPUT, --input INPUT
                        Input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        Base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
  -if {sdf,json}, --informat {sdf,json}
                        Input format. When using STDIN this must be specified.
  -of {sdf,json}, --outformat {sdf,json}
                        Output format. Defaults to 'sdf'.
  --meta                Write metadata and metrics files
```

```sh
zcat ../data/Kinase_inhibs.sdf.gz | python python/pipelines/rdkit/cluster_butina.py -t 0.6 -if sdf
```
