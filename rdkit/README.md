#RDKit examples

##1. screen.py
Performs ligand based virtual screening using RDKit descriptors

### Native execution

These examples need Python and RDKit to be present

```sh
python screen.py -h
usage: screen.py [-h] [-smiles SMILES] [-molfile MOLFILE] [-simmin SIMMIN]
                 [-simmax SIMMAX] [-d {maccs,morgan2,morgan3,rdkit}]
                 [-m {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}]
                 [-i INPUT] [-o OUTPUT]

RDKit screen

optional arguments:
  -h, --help            show this help message and exit
  -smiles SMILES        query structure as smiles (incompatible with -molfile
                        arg)
  -molfile MOLFILE      query structure as filename in molfile format
                        (incompatible with -smiles arg)
  -simmin SIMMIN        similarity lower cutoff (1.0 means identical)
  -simmax SIMMAX        similarity upper cutoff (1.0 means identical)
  -d {maccs,morgan2,morgan3,rdkit}, --descriptor {maccs,morgan2,morgan3,rdkit}
                        descriptor or fingerprint type (default rdkit)
  -m {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}, --metric {asymmetric,braunblanquet,cosine,dice,kulczynski,mcconnaughey,rogotgoldberg,russel,sokal,tanimoto}
                        similarity metric (default tanimoto)
  -i INPUT, --input INPUT
                        input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
```

```sh
zcat ../data/dhfr_3d.sdf.gz | python screen.py -smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2' -simmin 0.5 | gzip > screened_mols.sdf.gz
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed smiles. Gzip the resulting SDF and write to screened_mols.sdf.gz

```sh
zcat ../data/dhfr_3d.sdf.gz | python screen.py -molfile ../data/pyrimethamine.mol -simmin 0.5 -simmax 0.8 > screened_mols.sdf
```
Screen the structures in dhfr_3d.sdf.gz for similarity to the specifed file in molfile format of between 0.5 and 0.8. Write results to screened_mols.sdf

### Execution in Docker


##2 conformers.py
Generates 3D conformers

```sh
python conformers.py -h
usage: conformers.py [-h] [-n NUM] [-a ATTEMPTS] [-r RMSD] [-c {rmsd,tfd}]
                     [-t THRESHOLD] [-e EMIN] [-i INPUT] [-o OUTPUT]

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
                        input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        output is used as base name of the other files.
```

```sh
zcat ../data/Kinase_inhibs.sdf.gz | python conformers.py -n 5 > out.sdf
```
Tries to generate 5 conformers for each structure in the input file (may generate less if 2 conformers differ by less than the default RMSD) 

```sh
python conformers.py -i ../data/Kinase_inhibs.sdf.gz -n 5 -a 20 -c rmsd > out.sdf
```
Input specified using -i argument. Tries to generate 5 conformers using a maximum of 20 attempts, and clusters the conformers using RMSD distance.


##3. o3dAlign.py
Performs rigid 3D alignment of structures
```sh
python o3dAlign.py -h
usage: o3dAlign.py [-h] [-t THRESHOLD] [-n NUM] [-a ATTEMPTS] [-r RMSD]
                   [-e EMIN] [-i INPUT] [-o OUTPUT]
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
                        input SD file, if not defined the STDIN is used
  -o OUTPUT, --output OUTPUT
                        base name for output file (no extension). If not
                        defined then SDTOUT is used for the structures and
                        'o3dalign' is used as base name of the other files.
```


```sh
zcat ../data/dhfr_3d.sdf.gz | python o3dAlign.py ../data/pyrimethamine.mol -n 5 -t 10 > my_output.sdf
```
Aligns structures from STDIN to the specified molecule (3D), generating up to 5 conformers to align. Redirects output structures to my_output.sdf.

```sh
zcat ../data/dhfr_3d.sdf.gz | python o3dAlign.py ../data/pyrimethamine.mol -n 5 -t 10 -o screen1
```

```sh
python o3dAlign.py ../data/pyrimethamine.mol -i ../data/dhfr_3d.sdf.gz -o o3da
```
Aligns the 3D structures in the file dhfr_3d.sdf.gz to the 3D structure pyrimethamine.mol writing the output to o3da.sdf.gz

