#!/usr/bin/env bash

echo "Testing screen.py reading from STDIN and writing to STDOUT"
zcat ../data/dhfr_3d.sdf.gz | python python/pipelines/rdkit/screen.py\
  --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
  --simmin 0.5\
  --informat sdf > /dev/null || echo "FAILED"

echo "Testing conformers.py reading from STDIN and writing to STDOUT"
zcat ../data/Kinase_inhibs.sdf.gz | python python/pipelines/rdkit/conformers.py -n 2 -if sdf > /dev/null || echo "FAILED"

echo "Testing o3dAlign.py reading from STDIN and writing to STDOUT"
zcat ../data/Kinase_inhibs.sdf.gz | python python/pipelines/rdkit/o3dAlign.py\
  ../data/pyrimethamine.mol -n 2 -t 10 -if sdf > /dev/null || echo "FAILED"

echo "Testing filter.py reading from STDIN and writing to STDOUT"
zcat ../data/Kinase_inhibs.sdf.gz | python python/pipelines/rdkit/filter.py --hacmin 25 --hacmax 30 -if sdf > /dev/null || echo "FAILED"

echo "Testing cluster_butina.py reading from STDIN and writing to STDOUT"
zcat ../data/Kinase_inhibs.sdf.gz | python python/pipelines/rdkit/cluster_butina.py -t 0.6 -if sdf > /dev/null || echo "FAILED"

echo "Finished"

