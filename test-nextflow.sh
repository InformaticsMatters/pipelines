#!/bin/bash
# Simple manual tests for Nextflow workflows
# Before running make sure you have the latest images by running `./gradlew buildDockerImages`

set -e

echo 'Running screen+config'
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
cp src/nextflow/docking/smog.nsd.nf tmp/nextflow.nf
cp src/nextflow/docking/smog.nsd.config tmp/nextflow.config
gzip -c data/smog/DCP2_1.pdb > tmp/protein.pdb.gz
cp data/smog/confs.data.gz tmp/ligands.data.gz
cp data/smog/confs.metadata tmp/ligands.metadata
cd tmp
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow sh -c 'nextflow run nextflow.nf -c nextflow.config --score 100.0 -with-docker'
cd ..
sudo rm -rf tmp/*


echo 'Running PLI in squonk mode'
sudo rm -rf tmp/*
cp src/nextflow/docking/plip.nsd.nf tmp/nextflow.nf
cp src/nextflow/docking/plip.nsd.config tmp/nextflow.config
gzip -c data/smog/DCP2_1.pdb > tmp/protein.pdb.gz
cp data/smog/confs.data.gz tmp/ligands.data.gz
cp data/smog/confs.metadata tmp/ligands.metadata
cd tmp
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow sh -c 'nextflow run nextflow.nf -c nextflow.config --score 100.0 -with-docker'
cd ..
sudo rm -rf tmp/*


echo 'Running rDock in squonk mode'
sudo rm -rf tmp/*
cp src/nextflow/docking/rdock.nsd.nf tmp/nextflow.nf
cp src/nextflow/docking/rdock.nsd.config tmp/nextflow.config
cp data/hivpr.config.zip tmp/config.zip
cp data/dhfr_3d.data.gz tmp/ligands.data.gz
cp data/dhfr_3d.metadata tmp/ligands.metadata
cd tmp
docker run -it --rm -v $PWD:$PWD:z -w $PWD -v /var/run/docker.sock:/var/run/docker.sock informaticsmatters/nextflow sh -c 'nextflow run nextflow.nf -c nextflow.config --num_dockings 1 --limit 40 --chunk 5 -with-docker'
cd ..
sudo rm -rf tmp/*