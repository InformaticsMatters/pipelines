FROM busybox
MAINTAINER Tim Dudgeon <tdudgeon@informaticsmatters.com>

# An image to populate a Service Descriptor desitination directory
# (SD_DST), which is normally mounted when the image is run,
# with built-in Service Descriptor files from a source directory (SD_SRC).

ENV SD_SRC /sd-src
ENV SD_DST /sd-dst

RUN mkdir -p ${SD_SRC}

# Copy all potential Service Descriptors into the image...
COPY python/ ${SD_SRC}/python/
COPY nextflow/ ${SD_SRC}/nextflow/

# Remove anything that doesn't look like a Service Descriptor...
RUN rm -f `find ${SD_SRC} -type f -not -name "*.json" -not -name "*.yml" -not -name "*.yaml"`

# On execution copy files from source to destination...
CMD cp -R ${SD_SRC}/* ${SD_DST}
