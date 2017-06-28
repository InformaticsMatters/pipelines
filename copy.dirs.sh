#!/bin/sh

./gradlew dist
rm -rf ../squonk/data/testfiles/docker-services/pipelines
cp -r build/dist/pipelines ../squonk/data/testfiles/docker-services/ 
echo "Files copied to ../squonk/data/testfiles/docker-services/"

rm -rf ../squonk/docker/deploy/data/docker-services/pipelines
cp -r build/dist/pipelines ../squonk/docker/deploy/data/docker-services/
echo "Files copied to ../squonk/docker/deploy/data/docker-services/"

