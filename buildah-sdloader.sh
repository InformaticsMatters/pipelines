#!/bin/bash -x

# Just echo/log the Buildah version we have...
buildah version

# An image to populate a Service Descriptor destination directory
# (SD_DST), which is normally mounted when the image is run,
# with built-in Service Descriptor files from a source directory (SD_SRC).

SD_SRC=/sd-src
SD_DST=/sd-dst

CTR=`buildah from busybox`

# Some environment variables
# (expected to be available inside the final container image)
buildah config --env SD_SRC=/sd-src $CTR
buildah config --env SD_DST=/sd-dst $CTR

# Copy all potential Service Descriptors into the image...
buildah copy $CTR src/python/ ${SD_SRC}/python/
buildah copy $CTR src/nextflow/ ${SD_SRC}/nextflow/

# Remove anything that doesn't look like a Service Descriptor...
buildah run $CTR -- rm -f `find ${SD_SRC} -type f -not -name "*.json" -not -name "*.yml" -not -name "*.yaml"`

# On execution copy files from source to destination...
buildah config --cmd "cp -R ${SD_SRC}/* ${SD_DST}" $CTR

# Apply some annotations
buildah config --author "Tim Dudgeon <tdudgeon@informaticsmatters.com>" $CTR

# Commit the image to a name...
buildah commit $CTR informaticsmatters/rdkit_pipelines-sdloader