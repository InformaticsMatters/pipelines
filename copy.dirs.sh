#!/bin/sh

./gradlew dist
rm -rf ../squonk/data/testfiles/docker-services/pipelines ../squonk/data/testfiles/docker-services/nextflow
cp -r build/dist/pipelines ../squonk/data/testfiles/docker-services/ 
cp -r build/dist/nextflow ../squonk/data/testfiles/docker-services/
echo "Files copied to ../squonk/data/testfiles/docker-services/"

rm -rf ../squonk/docker/deploy/data/docker-services/pipelines ../squonk/docker/deploy/data/docker-services/nextflow
cp -r build/dist/pipelines ../squonk/docker/deploy/data/docker-services/
cp -r build/dist/nextflow ../squonk/docker/deploy/data/docker-services/
echo "Files copied to ../squonk/docker/deploy/data/docker-services/"

