#!/usr/bin/env bash
# run locally with something like this:
# ./post-service-descriptors.sh http://localhost:8091/coreservices/rest/v1/services
# or
# docker run -it --rm -v $PWD:$PWD:Z -w $PWD --network deploy_squonk_back centos:7 ./post-service-descriptors.sh

set -e

POST=${1:-http://coreservices:8080/coreservices/rest/v1/services}
BASE_D='docker://github.com/InformaticsMatters/pipelines'
BASE_N='nextflow://github.com/InformaticsMatters/pipelines'
CT_DJ="application/x-squonk-service-descriptor-docker+json"
CT_DY="application/x-squonk-service-descriptor-docker+yaml"
CT_MM="multipart/mixed"


for d in 'src/python/pipelines/dmpk' 'src/python/pipelines/docking' 'src/python/pipelines/rdkit' 'src/python/pipelines/dimorphite'
do
    for file in $d/*.dsd.yml
    do
	    echo $file
	    curl -X POST \
         -T $file\
         -H "Content-Type: $CT_DY"\
         -H "Base-URL: $BASE_D"\
         $POST
         echo ""
    done
done

for d in 'src/nextflow/docking' 'src/nextflow/rdkit'
do
    for file in $d/*.nsd.yml
    do
	    basename=${file::-4}
	    echo $basename
	    curl -X POST \
         -F "nextflow.nsd.yml=@${basename}.yml;type=application/x-squonk-service-descriptor-nextflow+yaml;filename=nextflow.nsd.yml"\
         -F "nextflow.nf=@${basename}.nf;type=text/plain;filename=nextflow.nf"\
         -F "nextflow.config=@${basename}.config;type=text/plain;filename=nextflow.config"\
         -H "Content-Type: $CT_MM"\
         -H "Base-URL: $BASE_N"\
         $POST
         echo ""
    done
done
