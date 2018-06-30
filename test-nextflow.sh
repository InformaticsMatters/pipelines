#!/bin/bash
# Simple manual tests for Nextflow workflows
# Before running make sure you have the latest images by running `./gradlew buildDockerImages`

set -e

echo 'Running screen+conformers in basic mode'
nextflow run src/nextflow/rdkit/screen+conformers.nf -c  src/nextflow/rdkit/screen.config -with-docker


echo 'Running SMoG2016 in basic mode'
nextflow run src/nextflow/docking/smog.nf -c src/nextflow/docking/smog.config -with-docker --ligands data/smog/confs.sdf --protein data/smog/DCP2_1.pdb informaticsmatters/smog:latest

echo 'Running PLI in basic mode'
nextflow run src/nextflow/docking/plip.nf -c src/nextflow/docking/plip.config -with-docker --ligands data/smog/confs.sdf --protein data/smog/DCP2_1.pdb informaticsmatters/pli:latest

echo 'Running rDock in basic mode'
nextflow run src/nextflow/docking/rdock.nf -c src/nextflow/docking/rdock.config -with-docker\
  --ligands data/hivpr_ligprep_100.sdf.gz\
  --protein data/hivpr_rdock.mol2\
  --asfile data/hivpr_rdock.as\
  --prmfile data/hivpr_rdock.prm\
  --num_dockings 2

echo 'Running SMoG2016 in squonk mode'
sudo rm -rf tmp/*
cd tmp
ln ../src/nextflow/docking/smog.nsd.nf nextflow.nf
ln ../src/nextflow/docking/smog.nsd.config nextflow.config
gzip -c ../data/smog/DCP2_1.pdb > protein.pdb.gz
ln ../data/smog/confs.data.gz ligands.data.gz
ln ../data/smog/confs.metadata ligands.metadata
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow-docker:0.30.2 sh -c 'nextflow run nextflow.nf -c nextflow.config --score 100.0 -with-docker'
cd ..

echo 'Running PLI in squonk mode'
sudo rm -rf tmp/*
cd tmp
ln ../src/nextflow/docking/plip.nsd.nf nextflow.nf
ln ../src/nextflow/docking/plip.nsd.config nextflow.config
gzip -c ../data/smog/DCP2_1.pdb > protein.pdb.gz
ln ../data/smog/confs.data.gz ligands.data.gz
ln ../data/smog/confs.metadata ligands.metadata
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow-docker:0.30.2 sh -c 'nextflow run nextflow.nf -c nextflow.config --score 100.0 -with-docker'
cd ..

echo 'Running rDock in squonk mode'
sudo rm -rf tmp/*
cd tmp
ln ../src/nextflow/docking/rdock.nsd.nf nextflow.nf
ln ../src/nextflow/docking/rdock.nsd.config nextflow.config
ln ../data/hivpr.config.zip config.zip
ln ../data/dhfr_3d.data.gz ligands.data.gz
ln ../data/dhfr_3d.metadata ligands.metadata
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow-docker:0.30.2 sh -c 'nextflow run nextflow.nf -c nextflow.config --num_dockings 1 --limit 40 --chunk 5 -with-docker'
cd ..

echo 'Running screen in squonk mode'
sudo rm -rf tmp/*
cd tmp
ln ../src/nextflow/rdkit/screen-dataset.nsd.nf nextflow.nf
ln ../src/nextflow/rdkit/screen-dataset.nsd.config nextflow.config
ln ../data/dhfr_3d.data.gz input.data.gz
ln ../data/dhfr_3d.metadata ligands.metadata
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow-docker:0.30.2\
  sh -c 'nextflow run nextflow.nf -c nextflow.config -with-docker --chunk 100 --simmin 0.5 --qsmiles "OC(=O)C1=CC=C(NC2=NC3=C(CN=C(C4=CC(Cl)=CC=C34)C3=C(F)C=CC=C3F)C=N2)C=C1"'
cd ..

echo 'Running screen-multi in squonk mode'
sudo rm -rf tmp/*
cd tmp
ln ../src/nextflow/rdkit/screen-multi-dataset.nsd.nf nextflow.nf
ln ../src/nextflow/rdkit/screen-multi-dataset.nsd.config nextflow.config
ln ../data/dhfr_3d.data.gz target.data.gz
ln ../data/nci100.data.gz query.data.gz
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow-docker:0.30.2\
  sh -c 'nextflow run nextflow.nf -c nextflow.config -with-docker --chunk 100 --simmin 0.55'
cd ..

sudo rm -rf tmp/*