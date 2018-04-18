#!/bin/bash -x

set -e

# Just echo/log the Buildah & Podman versions we have...
buildah version
podman version

# An image to populate a Service Descriptor destination directory
# (SD_DST), which is normally mounted when the image is run,
# with built-in Service Descriptor files from a source directory (SD_SRC).

SD_SRC=/sd-src
SD_DST=/sd-dst
TGT_IMG=informaticsmatters/rdkit_pipelines-sdloader
TGT_TAG=latest
REGISTRY=172.30.23.200

# Start a container image...
CNTR=$(buildah from busybox)

# Some environment variables
# (expected to be available inside the final container image)
buildah config --env SD_SRC=/sd-src ${CNTR}
buildah config --env SD_DST=/sd-dst ${CNTR}

# Copy all potential Service Descriptors into the image...
buildah run ${CNTR} -- mkdir -p ${SD_SRC}
buildah copy ${CNTR} src/python/ ${SD_SRC}/python/
buildah copy ${CNTR} src/nextflow/ ${SD_SRC}/nextflow/

# Remove anything that doesn't look like a Service Descriptor...
buildah run ${CNTR} -- rm -f `find ${SD_SRC} -type f -not -name "*.json" -not -name "*.yml" -not -name "*.yaml"`

# On execution copy files from source to destination...
buildah config --cmd "cp -R ${SD_SRC}/* ${SD_DST}" ${CNTR}

# Apply some annotations
buildah config --author "Tim Dudgeon <tdudgeon@informaticsmatters.com>" ${CNTR}

# Commit the image to a name...
buildah commit ${CNTR} ${TGT_IMG}:${TGT_TAG}

# List local images...
buildah images

# Push...
podman login --tls-verify=false --username jenkins --password $(oc whoami -t) ${REGISTRY}:5000
buildah push --tls-verify=false ${TGT_IMG}:${TGT_TAG} docker://${REGISTRY}:5000/${TGT_IMG}:${TGT_TAG}
podman logout ${REGISTRY}:5000
