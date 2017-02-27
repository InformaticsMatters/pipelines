#!/bin/sh

./gradlew dist
rm -rf ../lac/data/testfiles/docker-services/pipelines
cp -r build/dist/pipelines ../lac/data/testfiles/docker-services/ 

rm -rf ../lac/docker/deploy/data/docker-services/pieplines
cp -r build/dist/pipelines ../lac/docker/deploy/data/docker-services/
