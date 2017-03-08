#!/bin/sh

./gradlew dist
rm -rf ../lac/data/testfiles/docker-services/pipelines
cp -r build/dist/pipelines ../lac/data/testfiles/docker-services/ 
echo "Files copied to ../lac/data/testfiles/docker-services/"

rm -rf ../lac/docker/deploy/data/docker-services/pipelines
cp -r build/dist/pipelines ../lac/docker/deploy/data/docker-services/
echo "Files copied to ../lac/docker/deploy/data/docker-services/"

