#!/usr/bin/env bash

# A simple deployment script for OpenShift.
# It is assumed that your OpenShift Squonk application has been deployed
# and that you are running this from the Squonk project/namespace
# and have logged in as an appropriate user.

set -e pipefail

oc create -f post-service-descriptors.yaml

echo "Squonk pipeline deployment is underway..."
